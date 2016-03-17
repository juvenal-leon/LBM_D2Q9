clear all
load /Users/juvenal/Documents/Tesis/LBM_D2Q9/LBM_D2Q9/ux.txt
load /Users/juvenal/Documents/Tesis/LBM_D2Q9/LBM_D2Q9/uy.txt
load /Users/juvenal/Documents/Tesis/LBM_D2Q9/LBM_D2Q9/solid.txt

figure
colormap(gray(2))
image(2-solid)
hold on
quiver(ux,uy)