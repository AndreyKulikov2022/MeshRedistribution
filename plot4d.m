function h_pl=plot4d(x_d4,z_d4)
X=squeeze(x_d4(1,:,:,:));
Y=squeeze(x_d4(2,:,:,:));
Z=round(squeeze(z_d4),13);
h_pl=surf(X,Y,Z);
end