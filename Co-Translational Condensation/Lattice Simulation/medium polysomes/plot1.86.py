import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, least_squares
import os
from matplotlib.lines import Line2D
from functools import partial
from numpy.linalg import inv


N_max =5000
N_nc=9 #number of side chains, 5,9,15 S M L

C_µM = 10**7/(6.02*27)  # Extra 10^6 is to convert M to µM
BS=200
V_B=BS**3

colors = ['tab:green', 'tab:olive', 'tab:purple', 'tab:brown', 'tab:grey', 'black']
#Sequence-related color coding
seq_colors = ['tab:cyan','tab:red','tab:brown', 'tab:orange','tab:green']
seq_names = ['Homo','N-term', 'C-term', 'M-Block', 'Cross']
seq_c = [0.027390,0.027296,0.027300,0.027215,0.027316]
R = [0.105,0.082,0.058,0.042,0.039,0.036]
# Define fitting function: y = c + a / (x^(1/3))
def radius(n,c): # in lattice unit
    return (3*n/c/4/np.pi)**(1/3)

def c_sat(x, a, c):
    return c + a / (x ** (1/3))


def c_sat2(x, a, b, c):
    return c + a / (x ** (1/3)) + b / (x ** (2/3))

# Define the fitting function
def binding_frac_curve(x, a, n_pro, sigma):
    b = n_pro ** (2 / 3) / sigma
    return (x + b + a - np.sqrt((x + b + a)**2 - 4*b*x)) / 2/x

def eps(a,N_nc,V):
    if N_nc==5:
        v=160
    elif N_nc==9:
        v=208
    elif N_nc==15:
        v=265
    return -np.log(V/v/a)+1

# --- Your model (vectorized global) ---
def binding_curve_global(x_and_npro, a, sigma):
    x, npro = x_and_npro
    b = (npro ** (2/3)) / sigma
    inside = (x + b + a)**2 - 4*b*x
    return (x + b + a - np.sqrt(inside)) / 2

# --- Convenience: residuals (absolute or relative) ---
def residuals_abs(theta, x, npro, y):
    a, sigma = theta
    return binding_curve_global((x, npro), a, sigma) - y

def residuals_rel(theta, x, npro, y, eps=1e-12):
    yhat = binding_curve_global((x, npro), *theta)
    return (yhat - y) / (np.abs(y) + eps)

# === 1) Global fit (choose abs or relative) ===
def fit_global(x, npro, y, mode="abs", p0=(1.0, 1.0)):
    if mode == "rel":
        res = least_squares(residuals_rel, x0=np.array(p0), args=(x, npro, y), max_nfev=20000)
        a_fit, sigma_fit = res.x
        # J ~ jacobian of residuals; covariance approx: s^2 * (J^T J)^(-1)
        dof = len(y) - 2
        s2 = np.sum(res.fun**2) / max(dof, 1)
        JTJ_inv = inv(res.jac.T @ res.jac)
        pcov = s2 * JTJ_inv
    else:
        popt, pcov = curve_fit(binding_curve_global, (x, npro), y, p0=p0, maxfev=20000)
        a_fit, sigma_fit = popt
    return (a_fit, sigma_fit, pcov)

# === 2) Residual plot ===
def plot_residuals(x, npro, y, a_fit, sigma_fit, title="Residuals"):
    yhat = binding_curve_global((x, npro), a_fit, sigma_fit)
    res = y - yhat
    fig, ax = plt.subplots(1,2, figsize=(8,3.5), constrained_layout=True)
    ax[0].scatter(yhat, res, s=12, alpha=0.7)
    ax[0].axhline(0, color='k', lw=1)
    ax[0].set_xlabel("Fitted value")
    ax[0].set_ylabel("Residual")
    ax[0].set_title("Residuals vs fit")

    ax[1].scatter(x, res, s=12, alpha=0.7)
    ax[1].axhline(0, color='k', lw=1)
    ax[1].set_xlabel("#mRNA (x)")
    ax[1].set_title("Residuals vs x")
    fig.suptitle(title)
    plt.savefig(f"{surface_binding}/seq_{i}_residuals.svg", format='svg', bbox_inches='tight', pad_inches=0.03)
    plt.close()


# === 3) SE/CI and correlation ===
def summarize_cov(pcov, param_names=("a","sigma")):
    se = np.sqrt(np.diag(pcov))
    corr = pcov / np.outer(se, se)
    for i, name in enumerate(param_names):
        print(f"{name}: SE = {se[i]:.4g}")
    print("Parameter correlation matrix:")
    print(np.round(corr, 3))
    return se, corr

# === 4) Profile likelihood for 'a' ===
def profile_likelihood_a(x, npro, y, a_hat, sigma_hat, grid_factor=20, span_log10=1.0, mode="abs"):
    """
    Build a grid for 'a' on a log scale around a_hat (±span_log10 decades),
    refit sigma for each fixed 'a', compute Δχ² and CI where Δχ² <= 3.84.
    """
    # compute RSS_min and sigma^2_hat
    if mode == "rel":
        # chi2 = sum(res_rel^2); already dimensionless
        res0 = residuals_rel((a_hat, sigma_hat), x, npro, y)
        rss_min = np.sum(res0**2)
        dof = len(y) - 2
        s2_hat = 1.0  # already normalized residuals
    else:
        yhat0 = binding_curve_global((x, npro), a_hat, sigma_hat)
        res0 = y - yhat0
        rss_min = np.sum(res0**2)
        dof = len(y) - 2
        s2_hat = rss_min / max(dof, 1)  # variance estimate

    # a grid (log around a_hat); guard against nonpositive a_hat
    a_center = max(a_hat, 1e-8)
    a_grid = np.logspace(np.log10(a_center) - span_log10,
                         np.log10(a_center) + span_log10, grid_factor)

    rss_list = []
    sigma_list = []

    # refit sigma with a fixed a
    for a_fixed in a_grid:
        def obj_sigma(sig):
            # enforce positivity via transform if needed
            theta = (a_fixed, sig[0])
            if mode == "rel":
                r = residuals_rel(theta, x, npro, y)
            else:
                r = residuals_abs(theta, x, npro, y)
            return r

        # initial guess near sigma_hat
        res = least_squares(obj_sigma, x0=np.array([sigma_hat]), bounds=(1e-12, np.inf), max_nfev=20000)
        sigma_refit = res.x[0]
        rss = np.sum(res.fun**2)
        sigma_list.append(sigma_refit)
        rss_list.append(rss)

    rss_arr = np.array(rss_list)
    delta_chi2 = (rss_arr - rss_min) / s2_hat

    # 95% CI where Δχ² <= 3.84 (1 dof)
    mask = delta_chi2 <= 3.84
    ci = (np.min(a_grid[mask]) if np.any(mask) else np.nan,
          np.max(a_grid[mask]) if np.any(mask) else np.nan)

    # plot profile
    plt.figure(figsize=(4.5,4))
    plt.semilogx(a_grid, delta_chi2, '-o', ms=4)
    plt.axhline(3.84, color='r', ls='--', label='95% cutoff (Δχ²=3.84)')
    plt.axvline(a_hat, color='k', ls=':', label=r'$\hat a$')
    if np.all(np.isfinite(ci)):
        plt.axvspan(ci[0], ci[1], color='gray', alpha=0.2, label='95% CI for a')
    plt.xlabel(r'$a=Ve^{\epsilon-1}/v$')
    plt.ylabel(r'Profile $\Delta\chi^2$')
    plt.legend()
    # plt.title('Profile likelihood for a')
    plt.tight_layout()
    plt.savefig(f"{surface_binding}/seq_{i}_profile_likelyhood.svg", format='svg', bbox_inches='tight', pad_inches=0.03)
    plt.close()
    return a_grid, np.array(sigma_list), delta_chi2, ci

# === 5) (Optional) Bootstrap ===
def bootstrap_params(x, npro, y, n_boot=300, mode="abs", p0=(1.0,1.0), random_state=0):
    rng = np.random.default_rng(random_state)
    n = len(y)
    boots = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)   # resample indices
        xb, nb, yb = x[idx], npro[idx], y[idx]
        try:
            a_b, s_b, _ = fit_global(xb, nb, yb, mode=mode, p0=p0)
            boots.append((a_b, s_b))
        except Exception:
            continue
    boots = np.array(boots)
    if len(boots) == 0:
        print("Bootstrap produced no successful fits.")
        return None
    # summarize
    a_ci = np.percentile(boots[:,0], [2.5, 50, 97.5])
    s_ci = np.percentile(boots[:,1], [2.5, 50, 97.5])
    print(f"Bootstrap a:   2.5%={a_ci[0]:.3g}, med={a_ci[1]:.3g}, 97.5%={a_ci[2]:.3g}")
    print(f"Bootstrap sigma: 2.5%={s_ci[0]:.3g}, med={s_ci[1]:.3g}, 97.5%={s_ci[2]:.3g}")

    # quick scatter
    plt.figure(figsize=(4.5,4))
    plt.scatter(boots[:,0], boots[:,1], s=10, alpha=0.4)
    plt.axvline(a_ci[1], color='k', ls=':')
    plt.axhline(s_ci[1], color='k', ls=':')
    plt.xlabel('a')
    plt.ylabel('sigma')
    plt.title('Bootstrap distribution')
    plt.tight_layout()
    plt.savefig(f"{surface_binding}/seq_{i}_Bootstrap.svg", format='svg', bbox_inches='tight', pad_inches=0.03)
    plt.close()
    return boots


# Collect results
c_sat_fit_results = []
rdf_fit_results = []
binding_fit_results = []



# Store results
surface_binding = f"Surface_Binding"
os.makedirs(surface_binding, exist_ok=True)

sat_conc = f"Saturation_Concentration"
os.makedirs(sat_conc, exist_ok=True)

pls_fraction = f"Polysome_Fraction"
os.makedirs(pls_fraction, exist_ok=True)

frac_vs_mRNA_1000, frac_vs_mRNA_N_max =[],[]

def plot(i):
    file_path = f"phase_seq_{i}.csv"

    df = pd.read_csv(file_path)

    N_column = "N_Protein"  # N pro in dense phase
    Ndil_column = "New_Protein"  # Defined as #Protein - N_Protein, in dilute phase
    N_pls_column = "N_mRNA" #pls means polysome
    dense_column = "Dense%"  # Column for Dense%
    Nt_pls_column = "#mRNA"  # First grouping column, initial #mRNA
    Nt_column = "#Protein"  # Second grouping column, initial #mPro
    frac_column = "N_mRNA/#mRNA"     # fraction of mRNA in droplet
    # Convert columns to numeric
    df[Nt_column] = pd.to_numeric(df[Nt_column], errors='coerce')
    df[N_column] = pd.to_numeric(df[N_column], errors='coerce')
    df[dense_column] = pd.to_numeric(df[dense_column], errors='coerce')
    df[Nt_pls_column] = pd.to_numeric(df[Nt_pls_column], errors='coerce')
    df[Nt_column] = pd.to_numeric(df[Nt_column], errors='coerce')
    # N_max =max(df[Nt_column])
    plt.figure(figsize=(8, 6))
    df[frac_column] = df[N_pls_column]  / df[Nt_pls_column]

    df_new = df[df[dense_column] > 0.2].dropna(subset=[N_column, frac_column])



    ### Surface binding
    # --- Step 1: Collect all data ---
    all_xdata = []
    all_ydata = []
    all_npro = []

    for _, group in df.groupby(Nt_column):  # Nt_column = "#Protein"
        all_xdata.extend(group[Nt_pls_column].values)     # x = #mRNA
        all_ydata.extend(group[N_pls_column].values)     # y = N_mRNA
        all_npro.extend(group[N_column].values)           # N_Protein

    all_xdata = np.array(all_xdata)
    all_ydata = np.array(all_ydata)
    all_npro = np.array(all_npro)

    ##old fit and residuals plot
    # try:
    #     popt, _ = curve_fit(binding_curve_global, (all_xdata, all_npro), all_ydata,
    #                 p0=[1.0, 1.0], maxfev=20000)
    #     a_fit, sigma_fit = popt
    #     eps_fit=eps(a_fit,N_nc,V_B)
    #     sigma_nm2=3*3*(4*np.pi)**(1/3)*(3/seq_c[i-1])**(2/3)*sigma_fit #in units of nm^2

    #     #How good is the fit? coefficient of determination
    #     residuals = all_ydata - binding_curve_global((all_xdata, all_npro), *popt)
    #     plt.scatter(all_xdata, residuals, alpha=0.6)
    #     plt.axhline(0, color='k', linestyle='--')
    #     plt.xlabel("xdata")
    #     plt.ylabel("Residuals")
    #     plt.title("Residuals of Fit")
    #     plt.savefig(f"{surface_binding}/seq_{i}_residuals.svg", format='svg', bbox_inches='tight', pad_inches=0.03)
    #     plt.close()

    #     ss_res = np.sum(residuals**2)
    #     ss_tot = np.sum((all_ydata - np.mean(all_ydata))**2)
    #     r_squared = 1 - (ss_res / ss_tot)
    #     print(f"Surface binding fit: R² = {r_squared:.4f}, a = {a_fit:.2f}, eps = {eps_fit:.2f}, sigma = {sigma_nm2:.0f}")
    #     binding_fit_results.append({"Seq": i,"R²": f"{r_squared:.2f}",  "a_opt": f"{a_fit:.2f}","eps_opt": f"{eps_fit:.2f}", "sigma_opt": f"{sigma_fit:2f}", "sigma_nm2": f"{sigma_nm2:0f}"})
    # except Exception as e:
    #     print(f"Global fit failed: {e}")
    #     raise


    # 1) Global fit (pick "abs" or "rel")
    a_fit, sigma_fit, pcov = fit_global(all_xdata, all_npro, all_ydata, mode="abs", p0=(1.0, 1.0))
    eps_fit=eps(a_fit,N_nc,V_B)
    sigma_nm2=3*3*(4*np.pi)**(1/3)*(3/seq_c[i-1])**(2/3)*sigma_fit #in units of nm^2
    print(f"Surface binding global fit: a={a_fit:.4g}, sigma={sigma_fit:.4g}, eps = {eps_fit:.2f}, sigma_nm2 = {sigma_nm2:.0f}")
    binding_fit_results.append({"Seq": i, "a_opt": f"{a_fit:.2f}","eps_opt": f"{eps_fit:.2f}", "sigma_opt": f"{sigma_fit:.2f}", "sigma_nm2": f"{sigma_nm2:.0f}"})

    # 2) Residuals
    plot_residuals(all_xdata, all_npro, all_ydata, a_fit, sigma_fit, title="Surface-binding model")

    # 3) SEs & correlation
    se, corr = summarize_cov(pcov, param_names=("a","sigma"))

    # 4) Profile likelihood for 'a' (Tyler’s #3)
    agrid, sigmas_profile, dchi2, ci95 = profile_likelihood_a(
        all_xdata, all_npro, all_ydata,
        a_hat=a_fit, sigma_hat=sigma_fit,
        grid_factor=25, span_log10=1.0, mode="abs"
    )
    print(f"Profile-likelihood 95% CI for a: {ci95}")

    # 5) (optional) Bootstrap
    _ = bootstrap_params(all_xdata, all_npro, all_ydata, n_boot=400, mode="abs", p0=(a_fit, sigma_fit))


    # Plot separate curves by N_Pro 
    for Nt in sorted(df[Nt_column].unique()):
        plt.figure(figsize=(4.5, 4))
        Nt_list = list(sorted(df[Nt_column].unique()))
        Nt_list = Nt_list[::-1]  # reverse for top-down order
        for idx, Nt in enumerate(Nt_list):
            group = df[df[Nt_column] == Nt]

            xdata = group[Nt_pls_column].values  # #mRNA
            ydata = group[N_pls_column].values  # N_mRNA
            n_pro = np.mean(group[N_column].values)

            # Plot data
            #plt.scatter(xdata, ydata, color=f"C{idx}", s=10, alpha=0.4)
            # Compute mean and std of ydata (N_mRNA in droplet) for each #mRNA
            grouped = group.groupby(Nt_pls_column).agg(
                mean_y=(N_pls_column, 'mean'),
                std_y=(N_pls_column, 'std'),
                count=(N_pls_column, 'count')
            ).reset_index()
            grouped['stderr_y'] = grouped['std_y'] / np.sqrt(grouped['count'])  # optional: std error

            # Plot uncertainty as error bars
            plt.errorbar(grouped[Nt_pls_column], grouped['mean_y'],
                         yerr=grouped['stderr_y'], fmt='o', color=colors[idx],
                        capsize=3, markersize=4, alpha=0.8,zorder=4)
            plt.errorbar(0, 0, 
                      yerr=0, color=colors[idx],
                     fmt='o-', markersize=4, capsize=3, alpha=0.7#, label=f"$N_{{pro}}$={int(Nt)}") 
                    , label=f"R={radius(n_pro,seq_c[i-1])*0.003:.3f} µm") 
            # legend_elements.append(legend[0])
            # Global fitted curve
            x_fit = np.linspace(min(xdata), max(xdata)+5, 200)
            b = n_pro ** (2 / 3) / sigma_fit
            y_fit = (x_fit + b + a_fit - np.sqrt((x_fit + b + a_fit)**2 - 4 * b * x_fit)) / 2

            plt.plot(x_fit, y_fit, color=colors[idx], linestyle='-')


        # plt.legend(loc='upper left', fontsize=14, frameon=False)
        plt.legend(
            title='Condensate size',
            fontsize=12,
            loc='upper left',
            bbox_to_anchor=(-0.03, 1.03),  # 往左(-x) 和 往上(+y) 调
            frameon=False
        )        
        plt.xlim(0,51)
        plt.ylim(0,50)
        plt.xlabel("Total polysome number $N_T$", fontsize=13)
        plt.ylabel("#Surface-bound polysomes $N_B$", fontsize=13)
        plt.text(0.53, 0.91, f' {seq_names[i-1]} Seq fit:\n  ε = {eps_fit:.2f} ${{k_BT}}$ \n  σ = {sigma_nm2:.0f} nm${{^2}}$', transform=plt.gcf().transFigure, 
         fontsize=13, ha='left', va='top')
        #plt.title(f'{seq_names[i-1]} Seq fit:\n  ε = {eps_fit:.2f} ${{k_BT}}$, σ = {sigma_fit:.2f}')
        # plt.tick_params(axis='both', labelsize=10) 
        plt.tight_layout()
        plt.savefig(f"{surface_binding}/seq_{i}.svg", format='svg',bbox_inches='tight', pad_inches=0.03)
        plt.close()

    for idx, (name, group) in enumerate(df_new.groupby(Nt_pls_column)):
    
        # 先把 N 转换成半径 (单位 µm)
        group = group.copy()
        group['R'] = radius(group[N_column], seq_c[i-1]) * 0.003

        # 按 #Protein 分组，统计半径的均值和误差
        grouped = group.groupby(Nt_column).agg(
            mean_r=('R', 'mean'),
            std_r=('R', 'std'),
            mean_frac=(frac_column, 'mean'),
            std_frac=(frac_column, 'std'),
            count=('R', 'count')
        ).reset_index()

        grouped['stderr_frac'] = grouped['std_frac'] / np.sqrt(grouped['count'])  # y 方向误差
        grouped['stderr_r'] = grouped['std_r'] / np.sqrt(grouped['count'])        # x 方向误差

        # 画图
        plt.errorbar(
            grouped['mean_r'], grouped['mean_frac'],
            xerr=grouped['stderr_r'], yerr=grouped['stderr_frac'],
            fmt='o-', color=colors[idx+1], capsize=3, alpha=0.7,
            label=f"#Polysomes={name}"
        )


    # Formatting
    # plt.xlim(0, )
    plt.ylim(0, 1)
    plt.xlim(0, 0.11)
    plt.tick_params(axis='both', labelsize=12) 
    plt.xlabel("Condensate size $R$", fontsize=15)
    plt.ylabel("Fraction of CTC  $N_B/N_T$", fontsize=15)
    legend = plt.legend(title="Number of poplysomes", fontsize=15)
    legend.get_title().set_fontsize(15)    # plt.grid(True)
    plt.savefig(f"{pls_fraction}/polysome%_seq_{i}.svg", format='svg', bbox_inches='tight', pad_inches=0.03)
    plt.close()


   ### c_sat fit
    df= df[df[dense_column] > 0.2]
    # Plot N- and N+
    # Define new Ndil_column as Nt - N_column
    df[Ndil_column] = df[Nt_column] - df[N_column]
    df_noPS=df[df[Nt_pls_column]<1]

    # **Step 1: Global Fit for `c_opt`**
    if N_nc==5  or N_nc==10:
        popt_global, _ = curve_fit(c_sat, df_noPS[N_column], df_noPS[Ndil_column], maxfev=12000)
        _, c_opt = popt_global  # Extract the  fitted `c_opt` from base line
        c_sat_fit_results.append({"Seq": i, "polysomes": "Global", "a_opt": None, "c_opt": f"{c_opt:2f}"})
    elif N_nc==15:
        popt_global, _ = curve_fit(c_sat2, df_noPS[N_column], df_noPS[Ndil_column], maxfev=12000)
        _,_, c_opt = popt_global 
        c_sat_fit_results.append({"Seq": i, "polysomes": "Global", "a_opt": None,"b_opt": None, "c_opt": f"{c_opt:2f}"})

    #print(f"Global fit result for seq {i}: c_opt = {c_opt:.4f}")

    # ** Per-mRNA Grouping & Fitting**
    plt.figure(figsize=(5.3, 4))
    for idx, (Nt_pls, group) in enumerate(df.groupby(Nt_pls_column)):  

        # 复制数据，避免修改原 df
        group = group.copy()

        # 把 N 转换为半径（单位 µm）
        group['R'] = radius(group[N_column], seq_c[i-1]) * 0.003

        # 按 #Protein 分组计算均值和标准差
        grouped = group.groupby(Nt_column).agg(
            mean_r=('R', 'mean'),
            mean_x=(N_column, 'mean'),
            std_r=('R', 'std'),
            mean_y=(Ndil_column, 'mean'),
            std_y=(Ndil_column, 'std'),
            count=('R', 'count')
        ).reset_index()

        grouped['stderr_y'] = grouped['std_y'] / np.sqrt(grouped['count'])
        grouped['stderr_r'] = grouped['std_r'] / np.sqrt(grouped['count'])

        # 画数据误差条
        plt.errorbar(
            grouped['mean_r'], grouped['mean_y'] * C_µM / V_B,
            xerr=grouped['stderr_r'], yerr=grouped['stderr_y'] * C_µM / V_B,
            color=colors[idx], fmt='o', markersize=3, capsize=3, alpha=0.7
        )

        # 画 legend 示例
        plt.errorbar(
            0, -1, xerr=0, yerr=0,
            color=colors[idx], fmt='o-', markersize=3, capsize=3, alpha=0.7,
            label=f"{Nt_pls}"
        )
        # : Fit with Fixed `c_opt`, Only Optimize `a`**
        if N_nc==5 or N_nc==10:
            def c_sat_fixed_c(x, a):
                return c_sat(x, a, c_opt)
        if N_nc==15:
            def c_sat_fixed_c(x, a, b):
                return c_sat2(x, a, b, c_opt)

        popt, _ = curve_fit(c_sat_fixed_c, grouped['mean_x'], grouped['mean_y'], maxfev=12000)
        a_opt = popt[0]
        if N_nc==15:
            b_opt = popt[1]

        N_dense = (a_opt / 3) ** (3 / 4)
        N_dilute = c_sat(N_dense, a_opt, c_opt)
        N_tot = N_dense + N_dilute

        

        if N_nc==5 or N_nc==10:
            c_sat_fit_results.append({"Seq": i, "polysomes": f"{Nt_pls:.2f}", "a_opt": f"{a_opt:.2f}", "c_opt": f"{c_opt:.1f}","N_dense": f"{N_dense:.1f}", "N_dilute": f"{N_dilute:.1f}", "N_tot": f"{N_tot:.1f}"})

        if N_nc==15:
            c_sat_fit_results.append({"Seq": i, "polysomes": f"{Nt_pls:.2f}", "a_opt": f"{a_opt:.2f}","b_opt": f"{b_opt:.2f}", "c_opt": f"{c_opt:.1f}","N_dense": f"{N_dense:.1f}", "N_dilute": f"{N_dilute:.1f}", "N_tot": f"{N_tot:.1f}"})

        if N_nc==5 or N_nc==10:
            #plt.scatter(radius(N_dense,seq_c[i-1])*0.003, N_dilute*C_µM/V_B, marker='*', color=colors[idx])
            # Generate x values for the fitted curve
            x_fit = np.linspace(110, 6*N_max, 1000)
            y_fit = c_sat(x_fit, a_opt, c_opt)
        if N_nc==15:
            # Generate x values for the fitted curve
            x_fit = np.linspace(110, 6*N_max, 1000)
            # x_fit = np.linspace(N_dense, 1.2*N_max, 200)
            y_fit = c_sat2(x_fit, a_opt, b_opt, c_opt)
        r_fit = radius(x_fit,seq_c[i-1])*0.003

        # Plot the fitted curve
        plt.plot(r_fit, y_fit*C_µM/V_B, linestyle="-", color=colors[idx])

    # Formatting
    plt.xlabel(f"Condensate size $R$ (µm)", fontsize=15)
    plt.ylabel("Dilute phase conc. $c^{Dilute}$ (µM)", fontsize=15)


    legend = plt.legend(title="#Polysomes $N_{T}$", fontsize=10,loc='lower left')
    legend.get_title().set_fontsize(13)
    # plt.xlim(0, 1.005*max(r_fit))
    plt.xlim(0, 0.182)
    plt.ylim(bottom=0,top=2.2)
    plt.yticks([0, 0.5, 1, 1.5, 2]) 
    plt.tick_params(axis='both', labelsize=10)
    plt.savefig(f"{sat_conc}/csat_r_seq_{i}.svg", format='svg', bbox_inches='tight', pad_inches=0.03)
    plt.close()



    plt.figure(figsize=(5.3, 4))
    for idx, (Nt, group) in enumerate(df.groupby(Nt_column)):  

        # 复制数据，避免修改原 df
        group = group.copy()

        # 把 N 转换为半径（单位 µm）
        group['R'] = radius(group[N_column], seq_c[i-1]) * 0.003

        # 按 #Polysome 分组计算均值和标准差
        grouped = group.groupby(Nt_pls_column).agg(
            mean_r=('R', 'mean'),
            mean_x=(N_column, 'mean'),
            std_r=('R', 'std'),
            mean_y=(Ndil_column, 'mean'),
            std_y=(Ndil_column, 'std'),
            count=('R', 'count')
        ).reset_index()

        grouped['stderr_y'] = grouped['std_y'] / np.sqrt(grouped['count'])
        grouped['stderr_r'] = grouped['std_r'] / np.sqrt(grouped['count'])

        # 画数据误差条
        plt.errorbar(
            grouped[Nt_pls_column], grouped['mean_y'] * C_µM / V_B,
            yerr=grouped['stderr_y'] * C_µM / V_B,
            color=colors[5-idx], fmt='o-', markersize=3, capsize=3, alpha=0.7
        )

        # 画 legend 示例
        plt.errorbar(
            0, -1, yerr=0,
            color=colors[5-idx], fmt='o-', markersize=3, capsize=3, alpha=0.7,
            label=f"R={R[5-idx]:.3f} µm"
        )


    # Formatting
    plt.xlabel("Total polysome number $N_T$", fontsize=15)
    plt.ylabel("Dilute phase conc. $c^{Dilute}$ (µM)", fontsize=15)


    legend = plt.legend(title="Condensate size $R$", fontsize=9,loc='lower left')
    # legend.get_title().set_fontsize(12)
    plt.xlim(0, 50)
    plt.ylim(bottom=0,top=2.2)
    plt.yticks([0, 0.5, 1, 1.5, 2]) 
    plt.tick_params(axis='both', labelsize=10)
    plt.savefig(f"{sat_conc}/csat_NT_seq_{i}.svg", format='svg', bbox_inches='tight', pad_inches=0.03)
    plt.close()

    def compute_frac_vs_mRNA(results,df, seq_id, N_cutoff, 
                         Nt_column="#Protein", Nt_pls_column="#mRNA"):

        df_cut = df[df[Nt_column] == N_cutoff]

        grouped = df_cut.groupby(Nt_pls_column).agg(
            mean_frac=(frac_column, 'mean'),
            std_frac=(frac_column, 'std'),
            count=(frac_column, 'count')
        ).reset_index()

        grouped['stderr_frac'] = grouped['std_frac'] / np.sqrt(grouped['count'])

        for _, row in grouped.iterrows():
            results.append({
                "Seq": seq_id,
                "polysomes": row[Nt_pls_column],
                "frac": row["mean_frac"],
                "stderr_frac": row["stderr_frac"]
            })


    compute_frac_vs_mRNA(frac_vs_mRNA_1000,df_new, seq_id=i, N_cutoff=1000)
    compute_frac_vs_mRNA(frac_vs_mRNA_N_max,df_new, seq_id=i, N_cutoff=N_max)

    
    file_path = f"droplet_seq_{i}.csv"

    rdf_df = pd.read_csv(file_path)

    N_column = "N_Protein"  # N pro in dense phase
    Size_column = "Size" #radius of drop in lattice unit
    # Size_column = Size_column*0.003 #radius of drop in micron

    p_opt_rdf, _ = curve_fit(radius, rdf_df[N_column], rdf_df[Size_column], maxfev=12000)
    c_opt_rdf = p_opt_rdf[0]
    rdf_fit_results.append({"Seq": i, "c_opt": f"{c_opt_rdf:2f}"})

    #plot the fitting function
    rdf_grouped = rdf_df.groupby('NPro').agg(
            mean_x=(N_column, 'mean'),
            std_x=(N_column, 'std'),
            mean_y=(Size_column, 'mean'),
            std_y=(Size_column, 'std'),
            count=(N_column, 'count')
        ).reset_index()

    rdf_grouped['stderr_y'] = rdf_grouped['std_y'] / np.sqrt(rdf_grouped['count'])

    # Plot data with error bars
    plt.errorbar(rdf_grouped['mean_x'], rdf_grouped['mean_y']*0.003, 
                 xerr=rdf_grouped['std_x'], yerr=rdf_grouped['stderr_y']*0.003, color=colors[idx],
                 fmt='o', markersize=3, capsize=3, alpha=0.7) 

    x_fit = np.linspace(0, N_max+100, 200)
    y_fit = radius(x_fit,c_opt_rdf)*0.003

    # Plot the fitted curve
    plt.xlabel('#Proteins in condensate $N_{Pro}$')
    plt.ylabel('Condensa size $R$')
    plt.xlim(0, N_max+100)
    plt.ylim(0, 0.119)
    plt.plot(x_fit, y_fit, linestyle="-", color=seq_colors[i-1],label=f"{seq_names[i-1]}")
    plt.legend()
    plt.savefig(f"N_r_{i}.svg", format='svg', bbox_inches='tight', pad_inches=0.03)
    plt.close()

# Run for all sequences
for i in range(1, 6):
    plot(i)



# Convert to DataFrame and save
df_frac = pd.DataFrame(frac_vs_mRNA_N_max)
df_frac.to_csv(f"frac_vs_polysome_{N_max}.csv", index=False)

df_frac = pd.DataFrame(frac_vs_mRNA_1000)
df_frac.to_csv("frac_vs_polysome_1000.csv", index=False)


# **Save fitting results**
binding_df = pd.DataFrame(binding_fit_results)
binding_df.to_csv("binding_fit_results.csv", index=False)

rdf_fit_df = pd.DataFrame(rdf_fit_results)
rdf_fit_df.to_csv("rdf_fit_results.csv", index=False)

results_df = pd.DataFrame(c_sat_fit_results)
results_df.to_csv("c_sat_fit_results.csv", index=False)



# plot mRNA fraction of certain condensate
def plot_frac_mRNA(N):
    # Load the data
    df = pd.read_csv(f"frac_vs_polysome_{N}.csv")

    # Convert columns to numeric
    df["Seq"] = pd.to_numeric(df["Seq"], errors='coerce')
    df["polysomes"] = pd.to_numeric(df["polysomes"], errors='coerce')
    df["frac"] = pd.to_numeric(df["frac"], errors='coerce')
    df["stderr_frac"] = pd.to_numeric(df["stderr_frac"],errors='coerce')
    binding_df["Seq"] = pd.to_numeric(binding_df["Seq"], errors='coerce')
    binding_df["a_opt"] = pd.to_numeric(binding_df["a_opt"], errors='coerce')
    binding_df["eps_opt"] = pd.to_numeric(binding_df["eps_opt"], errors='coerce')
    binding_df["sigma_opt"] = pd.to_numeric(binding_df["sigma_opt"], errors='coerce')
    count=0
    # Plot
    plt.figure(figsize=(6, 5))
    # Iterate over unique sequences
    #for seq in df["Seq"].unique():
    legend_elements=[]
    for seq in [1,4,2,3,5]:
        subset = df[df["Seq"] == seq]
        para_subset = binding_df[binding_df["Seq"] == seq]
        plt.errorbar(subset["polysomes"], subset["frac"],
                     yerr=subset['stderr_frac'], fmt='o', color=seq_colors[seq-1],
                    capsize=3, markersize=4, alpha=0.8,zorder=4)

        label_str = f'{seq_names[seq-1]:<6} $\\varepsilon$ = {para_subset["eps_opt"].values[0]:.2f} $k_BT$'
        plt.errorbar(-1, -1, 
                  yerr=0, color=seq_colors[seq-1],
                 fmt='o-', markersize=4, capsize=3, alpha=1, label=label_str) 

        N_list=np.linspace(0,51,1000)
        frac_fit=[binding_frac_curve(x, para_subset["a_opt"], N, para_subset["sigma_opt"]) for x in N_list]
        plt.plot(N_list,frac_fit, color=seq_colors[seq-1],linewidth=2)
        count+=1

    # Formatting
    plt.ylim(0,1)
    plt.xlim(0,51)
    plt.xlabel("Total polysome number $N_T$",fontsize=16)
    plt.ylabel(f"Fraction of CTC  $N_B/N_T$",fontsize=16)
    plt.tick_params(axis='both', labelsize=13)
    legend = plt.legend(
        fontsize=13,
        loc='lower left',
        bbox_to_anchor=(-0.01, -0.03),  # x=0 (left), y=-0.3 (further down)
        frameon=False
    )
    for text in legend.get_texts():
        text.set_fontfamily('monospace')
    legend.get_title().set_fontfamily('monospace')



    plt.savefig(f"frac_polysome_{N}.svg", format='svg', bbox_inches='tight', pad_inches=0.03)


plot_frac_mRNA(1000)
plot_frac_mRNA(5000)


def plot_frac_mRNA_bar(N):
    # Load the data
    df = pd.read_csv(f"frac_vs_polysome_{N}.csv")

    # Convert columns to numeric
    df["Seq"] = pd.to_numeric(df["Seq"], errors='coerce')
    df["polysomes"] = pd.to_numeric(df["polysomes"], errors='coerce')
    df["frac"] = pd.to_numeric(df["frac"], errors='coerce')
    df["stderr_frac"] = pd.to_numeric(df["stderr_frac"], errors='coerce')

    # Filter only polysomes = 10
    df_10 = df[df["polysomes"] == 10]

    # Keep order [1,4,2,3,5] same as in your line plot
    ordered_seqs = [1, 4, 2, 3, 5]
    df_10 = df_10.set_index("Seq").loc[ordered_seqs].reset_index()

    # Bar plot
    plt.figure(figsize=(6,5))
    bars = plt.bar(
        df_10["Seq"].astype(str),
        df_10["frac"],
        yerr=df_10["stderr_frac"],
        capsize=4,
        color=[seq_colors[s-1] for s in df_10["Seq"]],
        alpha=0.8,
        zorder=3
    )

    # Formatting
    plt.ylim(0, 1)
    plt.xlabel("Sequence", fontsize=16)
    plt.ylabel("Fraction of CTC $N_B/N_T$", fontsize=16)
    plt.tick_params(axis='both', labelsize=13)
    plt.grid(axis='y', linestyle='--', alpha=0.6, zorder=0)

    # Custom x tick labels
    plt.xticks(
        ticks=range(len(df_10)),
        labels=[seq_names[s-1] for s in df_10["Seq"]],
        rotation=0
    )

    plt.savefig(f"frac_polysome_{N}_bar_Nt10.svg", format='svg', bbox_inches='tight', pad_inches=0.03)
    plt.close()



plot_frac_mRNA_bar(1000)
plot_frac_mRNA_bar(5000)

# plot mRNA fraction of certain condensate
df = pd.read_csv("c_sat_fit_results.csv")

# Convert columns to numeric
df["Seq"] = pd.to_numeric(df["Seq"], errors='coerce')
df["polysomes"] = pd.to_numeric(df["polysomes"], errors='coerce')
df["N_dilute"] = pd.to_numeric(df["N_dilute"], errors='coerce')

# Plot
plt.figure(figsize=(8, 6))
count=0
# Iterate over unique sequences
for seq in df["Seq"].unique():
    subset = df[df["Seq"] == seq]
    plt.plot(subset["polysomes"], subset["N_dilute"], marker='o', label=f'Seq {int(seq)}',color=seq_colors[count])
    count+=1

# Formatting
plt.ylim(100,300)
plt.xlabel("Total polysome number $N_T$",fontsize=15)
plt.ylabel("#Proteins in dilute phase $N_c$",fontsize=15)
plt.legend(fontsize=12)

plt.savefig(f"N-_polysomes.svg", format='svg', bbox_inches='tight', pad_inches=0.03)

C_µM = 10**7/(6.02*27)  # Extra 10^6 is to convert M to µM
BS=200
V_B=BS**3


plt.figure(figsize=(8, 6))
count=0
# Iterate over unique sequences
for seq in df["Seq"].unique():
    subset = df[df["Seq"] == seq]
    plt.plot(subset["polysomes"], C_µM*subset["N_dilute"]/V_B, marker='o', label=f'Seq {int(seq)}',color=seq_colors[count])
    count+=1

# Formatting
plt.ylim(0, 2.5) 
plt.xlabel("Total polysome number $N_T$",fontsize=15)
plt.ylabel("Nucleation-limited saturation concentration $c_{sat}^*$ (µM)",fontsize=12)
plt.legend(fontsize=15)

plt.savefig(f"cSat_polysomes.svg", format='svg', bbox_inches='tight', pad_inches=0.03)