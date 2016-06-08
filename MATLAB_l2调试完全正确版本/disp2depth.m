function depth = disp2depth(disp,f,baseline)
   disp = disp/8;
   disp = ones(size(disp))*mean(disp(:));
   depth = (f*baseline*ones(size(disp)))./disp;
end