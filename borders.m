function [borders_x, borders_y] = borders(img)
c = 1 ;
jmax=size(img, 1)-1;
imax=size(img, 2)-1;
for j = 2 : jmax
    for i = 2 : imax
    if img(j, i)
      if ~img(j-1, i) || ~img(j+1,i) || ~img(j,i+1) || ~img(j,i-1)
        borders_x(c) = i ;
        borders_y(c) = j ;
        c = c + 1 ;
      endfor
    endfor
  endfor
endfor
      