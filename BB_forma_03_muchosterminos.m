1;
pkg load optim;
pkg load miscellaneous;
# Lectura de los puntos y reescritura con BB centrada
xx=load('bordegota');
xxmedio=mean(xx)
xr=xx(:,2)-xxmedio(:,2);
yr=xx(:,1)-xxmedio(:,1);

# Ahora complejo para sacar modulo y argumento
# más fácil

xxcomp=xr+i*yr;
mod_rel=abs(xxcomp);
arg_rel=arg(xxcomp);


# Aca transformo las variables en las que usa el fiteo
x=arg_rel;
ydata=mod_rel;

#función con la que se hará el fiteo
f = @(p, x) p(1)+p(2)*cos(x)+p(3)*cos(2*x)+p(4)*cos(3*x)+p(5)*cos(4*x)+p(6)*cos(5*x)+p(7)*cos(6*x)+p(8)*cos(7*x)+p(9)*cos(8*x)+p(10)*cos(9*x)+p(11)*cos(10*x)+p(12)*cos(11*x)+p(13)*cos(12*x)+p(14)*cos(13*x)+p(15)*cos(14*x)+p(16)*cos(15*x)+p(17)*cos(16*x)+p(18)*sin(x)+p(19)*sin(2*x)+p(20)*sin(3*x);
#valores iniciales para el fiteo
#init_cvg = [1; 1.1;1.2;1;1;1;1;1;1];
#ahora con polinomios de Legendre
#f = @(p, x) p(1)+p(2)*legendrepoly(1,cos(x))+p(3)*legendrepoly(2,cos(x))+p(4)*legendrepoly(3,cos(x))+p(5)*legendrepoly(4,cos(x))+p(6)*legendrepoly(5,cos(x))+p(7)*legendrepoly(6,cos(x))+p(8)*legendrepoly(7,cos(x))+p(9)*legendrepoly(8,cos(x))+p(10)*legendrepoly(9,cos(x))+p(11)*legendrepoly(10,cos(x))+p(12)*legendrepoly(11,cos(x))+p(13)*legendrepoly(12,cos(x))+p(14)*legendrepoly(13,cos(x))+p(15)*legendrepoly(14,cos(x))+p(16)*legendrepoly(15,cos(x))+p(17)*legendrepoly(16,cos(x))+p(18)*legendrepoly(17,cos(x))+p(19)*legendrepoly(18,cos(x))+p(20)*legendrepoly(19,cos(x))+p(21)*legendrepoly(20,cos(x))+p(22)*legendrepoly(21,cos(x))+p(23)*legendrepoly(22,cos(x))+p(24)*legendrepoly(23,cos(x))+p(25)*legendrepoly(24,cos(x))+p(26)*legendrepoly(25,cos(x))+p(27)*legendrepoly(26,cos(x))+p(28)*legendrepoly(27,cos(x))+p(29)*legendrepoly(28,cos(x))+p(30)*legendrepoly(29,cos(x));
#valores iniciales para el fiteo
init_cvg = ones(20 ,1) ;

#este es el fiteo
[pc, mod_valc, cvgc, outpc] = nonlin_curvefit(f, init_cvg, x, ydata);

# Muestro resultados de los parametros
disp('Valores resultantes: ')
disp(pc)
pc_transp=pc';
#save -append -ascii salidaBB pc_transp
save -append -ascii 30legendre_coeficientes pc_transp

% graficamos la burbuja

% borde a partir de la imagen
xxx=ydata.*cos(x);
yyy=ydata.*sin(x);

% borde a partir del ajuste
angulo_ajuste = min(x):0.01:max(x) ;
r_ajuste = f(pc, angulo_ajuste) ;
x_ajuste = r_ajuste .* cos(angulo_ajuste) ;
y_ajuste = r_ajuste .* sin(angulo_ajuste) ;

fig = figure(2)
plot(x_ajuste, y_ajuste, 'r', 'LineWidth', 3)
hold on
scatter(xxx, yyy, 'b', 'filled')
xlim([-90 90])
ylim([-90 90])
#axis equal
saveas (fig,'figura2.png')


ff = figure(1)
scatter(x, ydata, 'b', 'filled');
hold on;
plot(min(x):0.01:max(x), f(pc, min(x):0.01:max(x)),'r', 'LineWidth',3);
saveas (ff,'figura3.png')


#figure; plot(x, ydata - mod_valc, "b.");


