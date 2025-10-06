% MATLAB script: Polysome and then Protein

% Polysome
len_mRNA = 100; % Length of mRNA
len_protein = 18; % Length of the protein part

global bead_index branch_start
bead_index = 1; % Initial bead index
fileID = fopen('z_Connection.txt', 'w'); % Open file to write connections

% Write connections for mRNA
for i = 1:(len_mRNA - 1)
    fprintf(fileID, '%d,%d,1\n', i, i + 1);
end

% Update bead index
bead_index = len_mRNA + 1;
len_5UTR = 18; % Length of 5'UTR
interval = 3;
branch_index = (len_5UTR+interval):interval:(len_5UTR + 54); % Generate branch indices
branch_start = [];
% Write connections for branches and protein regions
for i = 1:length(branch_index)-3
    fprintf(fileID, '%d,%d,1\n', branch_index(i), bead_index); % Branch connection
    fprintf(fileID, '%d,%d,1\n', bead_index, bead_index + 1); % Initial branch protein connection
    bead_index = bead_index + 1;

    % Write connections for extended protein region
    branch_start = [branch_start,bead_index+2];
    for j = 0:(interval/3 * i)
        fprintf(fileID, '%d,%d,1\n', bead_index, bead_index + 1);
        bead_index = bead_index + 1;
    end
    bead_index = bead_index + 1; % Increment bead index for next region
end

% Write connections for the main protein
for i = 1:(len_protein - 1)
    fprintf(fileID, '%d,%d,1\n', bead_index, bead_index + 1);
    bead_index = bead_index + 1;
end

fclose(fileID); % Close the file

branch_start
% Seq1 homopolymer
interaction1=[0 0 0; 0 -350 0; 0 0 0]; 
protein_list1 = [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2];
eps_list1 = energy_list(protein_list1);
writematrix(eps_list1,'z_BeadType_1.txt');
writematrix(interaction1, 'z_Energy_1.txt');

% Seq2 N-polymer
interaction2=[0 0 0; 0 -2200 -400; 0 -400 -50];
protein_list2 = [2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3];
eps_list2 = energy_list(protein_list2);
writematrix(eps_list2,'z_BeadType_2.txt');
writematrix(interaction2, 'z_Energy_2.txt');

% Seq3 C-polymer
interaction3=[0 0 0; 0 -2200 -400; 0 -400 -50];
protein_list3 = [3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2];
eps_list3 = energy_list(protein_list3);
writematrix(eps_list3,'z_BeadType_3.txt');
writematrix(interaction3, 'z_Energy_3.txt');

% Seq4 Block-polymer
interaction4=[0 0 0; 0 -50 -100; 0 -100 -1000];
protein_list4 = [3,3,2,2,2,3,3,2,2,2,3,3,2,2,2,3,3,3];
eps_list4 = energy_list(protein_list4);
writematrix(eps_list4,'z_BeadType_4.txt');
writematrix(interaction4, 'z_Energy_4.txt');

% Seq5 Cross-polymer
interaction5=[0 0 0; 0 -50 -400; 0 -400 -50];
protein_list5 = [2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3];
eps_list5 = energy_list(protein_list5);
writematrix(eps_list5,'z_BeadType_5.txt');
writematrix(interaction5, 'z_Energy_5.txt');

% 
% % Seq6 Diblock2-polymer
% interaction6=[-200 -400; -400 -50]
% protein_list6 = [1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2];
% eps_list6 = energy_list(protein_list6);
% eps6 = energy_table('z_Energy_6.txt', eps_list6, interaction6)

length(branch_start)

function  energy_list = energy_list(protein_list)
    global bead_index branch_start
    energy_list = ones(bead_index, 1);
    for i = 1: length(branch_start)
        energy_list(branch_start(i):branch_start(i)+i-1) = flip(protein_list(1:i)); % Reverse the order of selection, because N terminal should be near the end
    end
    energy_list(bead_index-17:bead_index) = protein_list(1:18);
end

% 
% function eps = energy_table(filename, energy_list, interaction)
%     eps = zeros(length(energy_list), length(energy_list));
%     for i = 1:length(energy_list)
%         for j = 1:length(energy_list)
%             if energy_list(i) * energy_list(j) ~= 0
%                 eps(i, j) = interaction(energy_list(i), energy_list(j));
%             end
%         end
%     end
%     writematrix(eps, filename);
% end