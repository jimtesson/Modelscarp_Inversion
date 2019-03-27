function sr = scrock(h,lambda,beta,rho_rock)

% sr = scrock(h,lambda,beta,rho_rock)
%
%--------------------------------------------------------------------------
% Schlagenhauf A., Gaudemer Y., Benedetti L., Manighetti I., Palumbo L.,
% Schimmelpfennig I., Finkel R., Pou K.
% G.J.Int., 2010
%-------------------------- ? ---------------------------------------------
%
%----------------------------- scrock.m -----------------------------------
%
% attenuation in the direction of e (perpendicular to fault scarp of dip beta)
% as a function of:
%   h = position in direction e of the sample
%   Lambda = true attenuation length (g.cm-2) (for ex. 208 for neutrons),
%   beta = scarp dip (degrees)
%   rho_rock = density (g.cm-3) of the rock.
%--------------------------------------------------------------------------

h(h==0) = 1/1000 ; % prevents NaNs 

m = 2.3 ; %  Lal exponent
beta = beta*pi/180 ;

[theta,phi] = meshgrid(0:90,0:180) ;
theta = theta*pi/180 ;
phi = phi*pi/180 ;
dphi = pi/180 ;
dtheta = pi/180 ;

% Upslope part : phi = [0 pi] , theta = [B(phi) pi/2]

B = atan(tan(beta).*sin(phi)) ; % apparent dip of the scarp in direction phi

da = f(beta,theta,phi) ;
da = exp(-h*rho_rock*da/lambda) ;
da = da.*(sin(theta).^m).*cos(theta).*(theta > B) ;
da = da*dphi*dtheta ;
da = sum(da(:)) ;

sa = da *(m + 1)/(2*pi) ;

% Downslope part : phi = [pi 2*pi] , theta = [0 pi/2]

dv = f(beta,theta,phi+pi) ;
dv = exp(-h*rho_rock*dv/lambda) ;
dv = dv.*(sin(theta).^m).*cos(theta) ;
dv = dv*dphi*dtheta ;
dv = sum(dv(:)) ;

sv = dv*(m + 1)/(2*pi) ;

sr = sa + sv ;

end

function d = f(beta,theta,phi)
num = 1 ;
den = sin(theta).*cos(beta) - sin(beta).*cos(theta).*sin(phi) ;
d = abs(num./den) ;
end