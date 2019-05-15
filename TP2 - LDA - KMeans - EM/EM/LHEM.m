function [ ppik ] = LHEM(l,clase,media,var,pi)
ppik=0;
for i=1:l
    ppik=ppik+(mvnpdf(clase(i,:),media,var)*pi);
end
end

