function eg_illustration(rho,Qx)
rho = uint8((rho - min(rho,[],'all'))/(max(rho,[],'all')-min(rho,[],'all'))*255);
rho = ind2rgb(rho,parula(256));

Qx = repmat(Qx, [1 1 3]);  

mask = Qx == 0;

Qx(mask) = rho(mask);
imshow(Qx);

end