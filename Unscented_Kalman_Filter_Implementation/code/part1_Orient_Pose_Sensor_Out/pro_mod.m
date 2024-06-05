function f_xun = pro_mod(sigma_points, angVel, acc, dt,n)
    
    f_xun = [];
    
    for i=1:(2*n+1)
        [R,G] = rot_G_matrix(sigma_points(4:6,i));
        G_inv = pinv(G);
        g=[0;0;-9.81];
         
        x_dot = [sigma_points(7,i);sigma_points(8,i);sigma_points(9,i); 
               G_inv*(angVel-sigma_points(10:12,i)-sigma_points(16:18,i)); 
          g+R*(acc-sigma_points(13:15,i)-sigma_points(19:21,i)); 
          sigma_points(22:24,i); sigma_points(25:27,i)];

        f_xun(1:15,i) = sigma_points(1:15,i)+dt*x_dot;

    end
end