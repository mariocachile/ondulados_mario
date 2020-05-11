frames = 450 : 5 : 500 ;
for i = 1 : length(frames)
  frame = frames(i) ;
  img = imread(['/home/burbujas/Desktop/bb3thr/BB3-2311-thresholded-frame',num2str(frame),'.tif']) ;
  img_filled = imfill(img, 'holes') ;
  img_filled = bwareaopen(img_filled, 3000) ;
  [borders_x, borders_y] = borders(img_filled) ;
  dlmwrite(strcat('/home/burbujas/Desktop/bb3thr/BB3-2311-borders-frame',num2str(frame),'.csv'),[borders_x' borders_y'], ",")
end