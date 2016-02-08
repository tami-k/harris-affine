function [ r, c ] = harris_affine( im, thresh, sigma, radius, disp )

    im_height  = size(im,1);
    im_width   = size(im,2);

    sigma_0 = sigma;
    sigma_step  = 1.2;
    sigma_m    = 13;
    sigma_array = (sigma_step.^(0:sigma_m-1))*sigma_0;
    
    harris_pts = zeros(0,3);
    for i=1:sigma_m
    
     s_I = sigma_array(i);   % интеграц. масштаб
     s_D = 0.7*s_I;          % дифф масштаб, берем 70%
     
     
      % маска
        x  = -3*s_D:3*s_D;
        dx = x .* exp(-x.*x/(2*s_D*s_D)) ./ (s_D*s_D*s_D*sqrt(2*pi));
        dy = dx';

        % 
        Ix = conv2(im, dx, 'same');
        Iy = conv2(im, dy, 'same');

        % матрица автокоррел€ции
        g   = fspecial('gaussian',max(1,fix(6*s_I+1)), s_I);
        Ix2 = conv2(Ix.^2, g,  'same');
        Iy2 = conv2(Iy.^2, g,  'same');
        Ixy = conv2(Ix.*Iy, g, 'same');

        % ќтклик
        k = 0.07;
        res = (Ix2.*Iy2 - Ixy.^2) - k*(Ix2 + Iy2).^2;	

        % локальный максимум
        [r,c,max_local] = findLocalMaximum(res,3*s_I);%3*s_I

        %TODO: порог в зав-ти от максимума

        % мера выше порога
        [r,c] = find(max_local>=thresh);
       

      %%  harris_pts  %
    end
     
    
     pp = size(harris_pts);
     for k=1:pp
         x(0) = harris_pts(k, i);
         U(0) = zeros(3);
       
      % номализованный Ћапласиан
      laplace = zeros(im_height,im_width,sigma_m);
     for i=1:sigma_m
        s_L = sigma_array(i);   % scale
        laplace(:,:,i) = s_L*s_L*imfilter(img,fspecial('log', floor(6*s_L), s_L),'replicate');
     end
    %провер€ем, где у нас макс. Ћапласа
    n   = size(harris_pts,1);
    points = zeros(n,3);
     for i=1:n
        r = harris_pts(i,1);
        c = harris_pts(i,2);
        s = harris_pts(i,3);
        val = laplace(r,c,s);
        
     end
    end

