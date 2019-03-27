function Ss = scsurf(Z,H,Lambda,beta,gamma,rho_rock)

% Ss = scsurf(Z,H,Lambda,beta,gamma,rho_rock)
% 
%--------------------------------------------------------------------------
% Schlagenhauf A., Gaudemer Y., Benedetti L., Manighetti I., Palumbo L.,
% Schimmelpfennig I., Finkel R., Pou K.
% G.J.Int., 2010
%-------------------------- ? ---------------------------------------------
%
%--------------------------- scsurf.m -------------------------------------
%
% Calculates the scaling factor Ss for the exhumed samples
% as a function of:
%   Z = depth (cm) measured on the scarp (0 at surface, > 0 above),
%   H = height of the sarp (cm),
%   Lambda = the true attenuation length (g.cm-2) (for ex. 208 for neutrons), 
%   beta = scarp dip (degrees),
%   gamma = dip of upper eroded part of the scarp, above beta (degrees),
%   rho_rock = density (g.cm-3) of the rock.
%--------------------------------------------------------------------------

hz=H-Z; hz(hz==0) = 1/1000; % to avoid NaNs

m = 2.3 ; % Lal exponent
beta = beta*pi/180 ;
gamma = gamma*pi/180 ;

[theta,phi] = meshgrid(0:90,0:180) ;
theta = theta*pi/180 ;
phi = phi*pi/180 ;
dphi = pi/180 ;
dtheta = pi/180 ;

% Downslope aerial part : phi = [pi 2*pi] , theta = [0 pi/2]
dv = (sin(theta).^m).*cos(theta) ;
dv = dv*dphi*dtheta ;
dv = sum(dv(:)) ;


% Upslope aerial part : phi = [0 pi] , theta = [B(phi) pi/2]
B = atan(tan(beta).*sin(phi)) ; % apparent dip of the scarp in the direction phi

da = (sin(theta).^m).*cos(theta).*(theta > B) ;
da = da*dphi*dtheta ;
da = sum(da(:)) ;

% Rock : phi = [0 pi] , theta = [C(phi) B(phi)]
C = atan(tan(gamma).*sin(phi)) ; % apparent dip of upper part of the scarp in direction phi

dr = f(gamma,beta,theta,phi) ;
dr = exp(-(hz)*rho_rock*dr/Lambda) ;
dr = dr.*(sin(theta).^m).*cos(theta).*(theta < B).*(theta > C) ;
dr = dr*dphi*dtheta ;
dr = sum(dr(:)) ;

S_air = (da + dv)*(m + 1)/(2*pi) ;
S_rock = dr*(m + 1)/(2*pi) ;
Ss = S_air + S_rock ;

end

function d = f(alpha,beta,theta,phi)
if (beta-alpha==0), alpha=alpha-0.0001; end % to avoid sin(0)
num = sin(beta - alpha) ;
den = sin(theta).*cos(alpha) - sin(alpha).*cos(theta).*sin(phi) ;
d = abs(num./den) ;
end