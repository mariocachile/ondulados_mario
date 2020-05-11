dirdoctorado = '/Doctorado/' ;
dirdoctorado = 'I:\Lucas\' ;

%addpath([dirdoctorado, filesep, 'repo-doctorado']) 
addpath([dirdoctorado, filesep, 'repo-doctorado', filesep, 'scripts-matlab']) 
addpath([dirdoctorado, filesep, 'repo-doctorado', filesep, 'scripts-matlab', filesep, 'fullfig']) 

scale = 8 / 1622 ; % cm / pix

data = csvread([dirdoctorado, 'Ondulados', filesep, 'Brahim', filesep, 'number_wavelength_amplitudepktopk_boolanalyzeornot.csv']) ;
bubble_number = data(:, 1) ;
wall_wavelength_cm = data(:, 2) ;
amplitude_pktopk_cm = data(:, 3) ;
analyze_or_not = logical(data(:, 4)) ;

% only sinusoidal walls for the moment
analyze_or_not(bubble_number > 69 & bubble_number < 104) = false ;
analyze_or_not(bubble_number > 125) = false ;

n = length(analyze_or_not) ;

[mean_distance_pix, max_distance_pix, min_distance_pix, wavelength_pix, ...
    mean_top_velocity_pixframe, area_pixsqr, fps] = initialize_variables(NaN(n, 1)) ;
folder_save = [dirdoctorado, 'Ondulados', filesep, 'results', ...
    filesep, 'average-parameters-and-shape'] ;

for bbl = 1 : 158
    if analyze_or_not(bbl)

        folder = [dirdoctorado, 'Ondulados', filesep, 'Brahim', filesep, 'bulle' , num2str(bbl), filesep] ;
        timestamps = get_timestamps(folder, 'tif') ;

        bg = rot90(max(imread([folder, num2str(timestamps(1)), '.tif']), imread([folder, num2str(timestamps(end)), '.tif']))) ;
        if bbl == 13
            [left_wall, right_wall] = walls_positions(bg, 0.2) ;
        else
            [left_wall, right_wall] = walls_positions(bg, 0.3) ;
        end
        distance_walls = right_wall - left_wall ;

        wavelength_cm = wall_wavelength_cm(bbl) ;
        amplitude_pk_cm = amplitude_pktopk_cm(bbl) ;
        
        % fitting the positions of the walls by a sine function
        [beta_left, rsqr_left, modelfun] = fit_walls(left_wall, ...
            [nanmean(left_wall), 0, amplitude_pk_cm / scale / 2, ...
            2 * pi * scale / wavelength_cm, 0]', 'sin plus linear') ;
        [beta_right, rsqr_right] = fit_walls(right_wall, ...
            [nanmean(right_wall), 0, amplitude_pk_cm / scale / 2, ...
            2 * pi * scale / wavelength_cm, 0]', 'sin plus linear') ;
        [beta_distance, rsqr_distance] = fit_walls(distance_walls, ...
            [nanmean(distance_walls), 0, amplitude_pk_cm / scale, ...
            2 * pi * scale / wavelength_cm, 0]', 'sin plus linear') ;
        mean_distance_pix(bbl) = beta_distance(1) ;
        max_distance_pix(bbl) = beta_distance(1) + abs(beta_distance(3)) ;
        min_distance_pix(bbl) = beta_distance(1) - abs(beta_distance(3)) ;
        wavelength_pix(bbl) = beta_distance(4) ;

        if rsqr_left < 0.97 || rsqr_right < 0.97 || rsqr_distance < 0.97
            warning('Look carefully the result of the fit of the walls.')
        end

        [n_frame_corrected, truedt_sec] = time_correction_ids(timestamps, [20 7]) ;
        t_sec = n_frame_corrected * truedt_sec ;
        fps(bbl) = 1 / truedt_sec ;
        
        [x_cm_pix, y_cm_pix, y_top_pix, x_top_pix, y_top_subpix, area, ...
            y_bottom_subpix] = initialize_variables(NaN(length(timestamps), 1)) ;
        [Xb, Yb] = initialize_variables(cell(length(timestamps), 1)) ;
        is_ordered = false(length(timestamps), 1) ;
        
        bad_equals = [] ; % will contain the frame indices that should be discarded
        % due to acquisition errors (two consecutive frame equal or equal by parts)
        mkdir([folder, filesep, 'test']) ;

        for i = 1 : length(timestamps)
            img_original = rot90(imread([folder, num2str(timestamps(i)), '.tif'])) ;
            if i > 1
                if isequal(img_original, img_prev)%isequal_byparts(img, img_prev, 10)
                    warning(['Bubble ', num2str(bbl), ', images ',num2str(i-1),' and ',num2str(i),' are equal. Acquisition error.'])
                    bad_equals(end + 1) = i - 1 ; % this should not arrive frequently (it's for that reason that I don't preallocate)
                    writelog([dirdoctorado, 'Ondulados', filesep, 'Brahim', filesep, 'debug.txt'], ['Bubble ', num2str(bbl), ', images ', num2str(timestamps(i-1)), ' and ', num2str(timestamps(i)), ' are equal. Acquisition error.'])
                end
            end    
            img_prev = img_original ;
            if bbl == 13
                img = image_preproc(img_original, bg, 0.2, 7000, 'clear border', 'fill holes') ;
            else
                img = image_preproc(img_original, bg, 0.3, 7000, 'clear border', 'fill holes') ;
            end
            %img(bwlabel(img,4) ~= 1) = 0 ; % leaves only 4-connected
            %objects (already done in image_preproc when calling
            %bwopenarea)
            %img = bwareaopen() ;

            A = regionprops(img, 'Area', 'Centroid',  'PixelList') ;
            if not(isempty(A))
                [area(i), ind] = max([A.Area]) ; % analyzing the object with the maximum surface
                pixellist = A(ind).PixelList ;
                x_cm_pix(i) = A(ind).Centroid(1) ;
                y_cm_pix(i) = A(ind).Centroid(2) ;
                [x_top_pix(i), y_top_pix(i)] = top_coordinates(pixellist) ;
                [~, t0] = crossing(double(img_prev(:, round(x_top_pix(i)))), 1 : size(img_prev, 1), 120) ;
                if length(t0) == 4
                    y_top_subpix(i) = mean(t0(1:2)) ;
                    y_bottom_subpix(i) = mean(t0(3:4)) ;
                end
                [xb, yb, niter] = shapecontour(pixellist) ;
                if niter == size(pixellist, 1)
                    writelog([dirdoctorado, 'Ondulados', filesep, 'Brahim', filesep, 'debug.txt'], ['Bubble ', num2str(bbl), ' frame ', num2str(timestamps(i)), ' had to run shapecontour_notordered.'])
                    [xb, yb] = shapecontour_notordered(pixellist) ;
                else
                    is_ordered(i) = true ;
                end
                Xb{i} = xb ;
                Yb{i} = yb ;
                [min_dist_leftwall, indices_mindist_left] = min_distance_to_wall([xb, yb], left_wall, (1:length(left_wall))) ;
                [min_dist_rightwall, indices_mindist_right] = min_distance_to_wall([xb, yb], right_wall, (1:length(right_wall))) ;
            else
                [pixellist, xb, yb] = initialize_variables([]) ;
                [min_dist_leftwall, min_dist_rightwall, indices_mindist_left, indices_mindist_right, ind] = initialize_variables(NaN) ;
            end
            %disp([i, length(timestamps)])
            
            %{
            f = figure('visible', 'off') ;
            imshow(img_original)
            hold on
            plot(left_wall, 1:length(left_wall), 'r')
            plot(right_wall, 1:length(right_wall), 'g')
            plot(modelfun(beta_left, 1:length(left_wall)), 1:length(left_wall), 'b--')
            plot(modelfun(beta_right, 1:length(right_wall)), 1:length(right_wall), 'b--')
            if not(isempty(A))
                scatter(x_top_pix(i), y_top_pix(i), 'r', 'filled')
                scatter(x_top_pix(i), y_top_subpix(i), 'r')
                scatter(x_top_pix(i), y_bottom_subpix(i), 'w')
                scatter(x_cm_pix(i), y_cm_pix(i), 'g', 'filled')
                plot(xb, yb, 'g')
                scatter(xb(indices_mindist_left), yb(indices_mindist_left), 'b', 'filled')
                scatter(xb(indices_mindist_right), yb(indices_mindist_right), 'b', 'filled')
            end
            %fullfig(f)
            save_figure(f, [folder, filesep, 'test', filesep, num2str(timestamps(i)), '.png']) ;
            close(f)
            %}
            
        end
        
        [x_cm_pix, y_cm_pix, x_top_pix, y_top_pix, area, n_frame_corrected, ...
            t_sec, y_top_subpix] = remove_elements(bad_equals, x_cm_pix, ...
            y_cm_pix, x_top_pix, y_top_pix, area, n_frame_corrected, t_sec, y_top_subpix) ;
        params_topfit = fit_linparams(n_frame_corrected, y_top_pix, y_top_pix*0.1, {@(x)1 @(x)x}) ;
        mean_top_velocity_pixframe(bbl) = params_topfit(2) ;
        ordinate = params_topfit(1) ;
        area_pixsqr(bbl) = nanmean(area) ;
        
        %{
        f = figure('visible', 'off') ;
        plot(n_frame_corrected, y_top_pix)
        xlabel('Frame number (corrected by frame loss)')
        ylabel('Top position (pix)')
        set(gca, 'FontSize', 20)
        fullfig(f) ;
        mkdir([folder, filesep, 'plots']) ;
        saveas(f, [folder, filesep, 'plots', filesep,'top_position_vs_frame.tif']) ;
        close(f)
        
        f = figure('visible', 'off') ;
        plot(n_frame_corrected, area)
        xlabel('Frame number (corrected by frame loss)')
        ylabel('Area (pix sqr)')
        set(gca, 'FontSize', 20)
        fullfig(f) ;
        saveas(f, [folder, filesep, 'plots', filesep,'area_vs_frame.tif']) ;
        close(f)
        %}
        %mkdir([folder, filesep, 'results']) ;
        save([folder_save, filesep, 'bubble',num2str(bbl),'.mat'], ...
            'bbl', 'timestamps', 'left_wall', 'right_wall', 'distance_walls', ...
            'wavelength_cm', 'amplitude_pk_cm', 'beta_left', 'rsqr_left', ...
            'beta_right', 'rsqr_right', 'beta_distance', 'rsqr_distance', ...
            'n_frame_corrected', 'truedt_sec', 't_sec', 'bad_equals', ...
            'area', 'x_cm_pix', 'y_cm_pix', 'x_top_pix', 'y_top_pix', ...
            'y_top_subpix', 'y_bottom_subpix', 'params_topfit', 'Xb', 'Yb', ...
            'is_ordered') ;
        
        %save([folder, filesep, 'results', filesep,'bubble', num2str(bbl), ...
        %    '_pixel_coordinates_borders.mat'], 'bbl', 'timestamps', 'Xb', ...
        %    'Yb') ;
        
        disp(['Bubble ', num2str(bbl), ' analyzed.'])
    end
end

%scatter(sqrt(4*area_pixsqr * scale^2/pi), -mean_top_velocity_pixframe * scale .* fps', 40, mean_distance_pix, 'filled')

 
