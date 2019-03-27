function sd = scdepth(Z,H,Lambda,alpha,beta,gamma,rho_rock,rho_coll)

% sd = scdepth(Z,H,Lambda,alpha,beta,gamma,rho_rock,rho_coll)
%
%--------------------------------------------------------------------------
% Schlagenhauf A., Gaudemer Y., Benedetti L., Manighetti I., Palumbo L.,
% Schimmelpfennig I., Finkel R., Pou K.
% G.J.Int., 2010
%-------------------------- ? ---------------------------------------------
%
%--------------------------- scdepth.m ------------------------------------
%
% Calculates the scaling factor sd for the buried samples
% as a function of:
%   Z = depth (cm) measured on the scarp (0 at surface, < 0 underneath),
%   H = height of the scarp (cm),
%   Lambda = the true attenuation length (g.cm-2) (for ex. 208 for neutrons),
%   alpha = colluvium dip (degrees) ;
%   beta = scarp dip (degrees),
%   gamma = dip of upper eroded part of the scarp, above beta (degrees),
%   rho_rock = density (g.cm-3) of the rock,
%   rho_coll = density (g.cm-3) of the colluvium.
%--------------------------------------------------------------------------

Z(Z==0) = -1/1000 ; % prevents NaNs

m = 2.3 ; %  Lal exponent
alpha = alpha*pi/180 ;
beta = beta*pi/180 ;
gamma = gamma*pi/180 ;

[theta,phi] = meshgrid(0:90,0:180) ;
theta = theta*pi/180 ;
phi = phi*pi/180 ;
dphi = pi/180 ;
dtheta = pi/180 ;

% Downslope part of colluvium : phi = [pi 2*pi] , theta = [0 pi/2]

dv = f(alpha,beta,theta,phi+pi) ;
dv = exp(Z*rho_coll*dv/Lambda) ; % Z negative under the colluvium
dv = dv.*(sin(theta).^m).*cos(theta) ;
dv = dv*dphi*dtheta ;
dv = sum(dv(:)) ;

% Upslope part of colluvium : phi = [0 pi] , theta = [B(phi) pi/2]

B = atan(tan(beta).*sin(phi)) ; %  apparent dip of scarp in direction phi

da = f(alpha,beta,theta,phi) ;
da = exp(Z*rho_coll*da/Lambda) ; % Z negative under the colluvium
da = da.*(sin(theta).^m).*cos(theta).*(theta > B) ;
da = da*dphi*dtheta ;
da = sum(da(:)) ;

% Rock : phi = [0 pi] , theta = [C(phi) B(phi)]

C = atan(tan(gamma).*sin(phi)) ; % apparent dip of upper part of scarp in direction phi

dr = f(alpha,beta,theta,phi) ;
dr = exp(-(H - Z)*rho_rock*dr/Lambda) ; % H - Z > H
dr = dr.*(sin(theta).^m).*cos(theta).*(theta < B).*(theta > C) ;
dr = dr*dphi*dtheta ;
dr = sum(dr(:)) ;

sc = (da + dv)*(m + 1)/(2*pi) ;
sr = dr*(m + 1)/(2*pi) ;
sd = sc + sr ;

end

function d = f(alpha,beta,theta,phi)
num = sin(beta - alpha) ;
den = sin(theta).*cos(alpha) - sin(alpha).*cos(theta).*sin(phi) ;
d = abs(num./den) ;
end