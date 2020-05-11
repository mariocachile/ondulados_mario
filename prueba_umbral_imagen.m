 img_original = rot90(imread("burbuja01.tif")) ;
 %img_originalsinrotar = imread(burbuja01, '.tif']) ;
 img_prev = img_original ;
 %           if bbl == 13
  %              img = image_preproc(img_original, bg, 0.2, 7000, 'clear border', 'fill holes') ;
   %         else
    %            img = image_preproc(img_original, bg, 0.3, 7000, 'clear border', 'fill holes') ;
     %       end
     %img = image_preproc(img_original, bg, 0.2, 7000, 'clear border', 'fill holes') ;
     
    imshow(img_original);
    [xx,yy] = borders(img_original)
    