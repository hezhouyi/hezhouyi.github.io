#!/bin/bash
source_code='G_6.11_DNA_nuc.f90'
python_analysis='analysis_mRNA2.6.py'
mkdir -p Runkey

N1=(0 1)
N2=(3000)



for serial in {1..20}; #replica
do
for i in {1..5}; #sequence
#for i in 1;
do 

for i1 in "${N1[@]}" ; do
for i2 in "${N2[@]}"; do

# Calculate Ti, so that simulation time scales with system size
Ti=$(echo "1500 * $i1 + 100 * $i2" | bc -l)
max_Ti=1000000

# Ensure Ti does not exceed max_Ti
if [ "$(echo "$Ti > $max_Ti" | bc -l)" -eq 1 ]; then
    Ti=$max_Ti
fi

# Convert Ti to an integer
Ti_int=$(printf "%.0f" "$Ti")

# Calculate R = (i2/2)^(1/3)
R=$(echo "2.5 * e( l($i2) / 3 )" | bc -l)

# Round up R to the nearest integer (ceiling)
R_int=$(echo "$R + 0.999999" | bc | cut -d. -f1)

# Print results for verification
# echo "Ti: $Ti_int"
# echo "R: $R_int"

index="seq${i}_${i1}_${i2}_${serial}"
    cat > "./Runkey/y_runkey_${index}.txt" <<EOF
Tmax
300000
Ti
${Ti_int}
Trelax
0
BoxSize
400 400 400
NkN
2
Nk
$i1 $i2
MaxSumNk
100 2000
Vk
208 18
Vak
81 18
HasDNA
F
UseGrand
F
UseESelf
T
LSteric
10
ESelf
1000
ChemPot
1000
BeadTypeFile
z_BeadType.txt
EnergyFile
z_Energy.txt
ConnectionFile
z_Connection.txt
UseSpringLinkerEq
F
E_SpringLinkerVector
0
SpringLinkerEqFile
z_SpringLinkerEq.txt
StartImported
F
StartDroplet
F
StartDropletR
$R_int
StartDropletNk
0 $i2
StartSlab
F
StartSlabH
30
StartSlabNk
1 1000
WriteLatticeEnd
F
iMoveRot
100
iMoveTransI
100
iMoveTransII
100
iMoveSlitherI
0
iMoveSlitherII
100
iMoveClusterI
10
iMoveClusterII
0
iMoveGrand
0
E_PrintStep
1000
E_Sanity
10000
E_Energy
1000
E_BoundSites
0
E_BondTypes
0
E_Lattice
0
WriteBoundTo
F
E_Network
0
E_SelfLoop
0
E_ClusterHist
1000
E_LargestCluster
50
WriteCoreCluster
F
E_XYZDist_All
0
E_XYZDist_None
0
E_RadDist_All
0
E_RadDist_Cluster
100
E_XYZDist_Cluster
0
WriteIndiDist_Cluster
F
E_RG_All
0
E_RG_Cluster
0
Seed
0
END
EOF
	cat > phase_${index}.slurm <<EOF
#! /bin/sh
#SBATCH -J ${index}
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --mem-per-cpu=512 # memory in MB used by the job
#SBATCH --time=168:00:00 # max 168 hours, default 8 hours
#SBATCH --error=${index}.err
#SBATCH --output=${index}.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zhouyi.he@mailbox.tu-dresden.de

hostname  # prints on which node the job is executed

mkdir -p phase_${index}
cd phase_${index}


cp \$SLURM_SUBMIT_DIR/${source_code} . #copy the executable
cp \$SLURM_SUBMIT_DIR/Runkey/y_runkey_${index}.txt ./y_runkey.txt #copy files which are read in from the program
cp \$SLURM_SUBMIT_DIR/z_Energy_${i}.txt ./z_Energy.txt
cp \$SLURM_SUBMIT_DIR/z_BeadType_${i}.txt ./z_BeadType.txt
cp \$SLURM_SUBMIT_DIR/z_Connection.txt .
#cp \$SLURM_SUBMIT_DIR/z_SpringLinkerEq.txt .
cp \$SLURM_SUBMIT_DIR/${python_analysis} . #copy python analysis program


gfortran ${source_code} -O3
./a.out y_runkey.txt # execute the simulation
module load release/24.10
module load Anaconda3/2024.02-1
#python ${python_analysis}  #run analysis with python

EOF
    sbatch "phase_${index}.slurm"  # Submit the SLURM job
done
done
done
done
