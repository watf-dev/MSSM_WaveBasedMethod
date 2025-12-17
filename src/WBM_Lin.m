%% ======================= STEP I. PREPARE ===============================
clear classes; clear all; clc;

%% ====================== STEP II. INPUT DATA ============================

% ---------- II.1 Geometry: bounding box (for WBM wave sets) -------------
% Bounding box defines the rectangular domain used for the WBM wave sets.
Lx = 1.25;      
Ly = 1.0;       

% Initialize nodal coordinates for the four corners.
% The x-coordinates of the right side (nodes 2 and 3) will be overwritten
% by the B-spline end points computed below.
nodes = [ 0,  0;   
          1,  0;   % (temporary value)
          1,  1;   % (temporary value)
          0,  1];  

% ---------- II.2 Material / wave properties -----------------------------
elasticityModulus = 150000;      
density           = 1.225;       % kg/m^3
waveVelocity      = 339;         % m/s

frequency         = 800;         % Hz
omega             = 2*pi*frequency;

% --- small damping to avoid irregular frequencies -----------------------
% NOT strictly required
% Pluymers et al. (2007) note that introducing a small damping term 
% is a common remedy for avoiding irregular (near-singular) frequencies.  
% With η = 1e-6, the physical solution is essentially unchanged.
eta       = 1e-6; 
omega_eff = omega*(1 + 1i*eta); 

% ---------- II.3 Boundary conditions: all rigid -------------------------
% Deckers 2010 B-spline validation example uses fully rigid walls
BCs = ['v'; 'v'; 'v'; 'v'];

% ---------- II.4 Truncation + quadrature -------------------------------
% Truncation of the two wave sets (cf. standard 2-set WBM formulation).
n1  = 10;
n2  = 10;

nGP = 50;                  % # Gauss points per edge

% ---------- II.5 Point source ------------------------------------------
q0 = 1;                    % amplitude [Pa]
xQ = 0.25;
yQ = 0.25;

%  **************************************
%  ***  II.6 B-spline for Right Edge  ***
%  **************************************

% Degree and open knot vector of the cubic B-spline
B_degree      = 3;
B_knot_vector = [0, 0, 0, 0, 0.5, 1, 1, 1, 1];

% Control points of the B-spline (right boundary of the cavity)
B_control_points = [
    1.00, 1.00;
    1.16, 0.84;
    1.34, 0.50;
    1.16, 0.16;
    1.00, 0.00
];

% --- Evaluate spline to check extents -----------------------------------
u_check  = linspace(0, 1, 2001); % fine sampling for a smooth spline curve
x_curve  = zeros(size(u_check));
y_curve  = zeros(size(u_check));

for i = 1:numel(u_check)
    [C, ~] = evalBspline(u_check(i), B_control_points, B_degree, B_knot_vector);
    x_curve(i) = C(1);
    y_curve(i) = C(2);
end

% --- End points (u=0, u=1) are used to connect to straight edges --------
[C_top, ~] = evalBspline(0, B_control_points, B_degree, B_knot_vector);  
[C_bot, ~] = evalBspline(1, B_control_points, B_degree, B_knot_vector);  

x_spline_end = C_bot(1);     % ≈ 1.0 (right-most x of the physical cavity)

% Update the straight-edge part of the physical boundary:
% The right edge is represented by the B-spline curve between C_top & C_bot.
nodes = [ 0,             0;          % 1: bottom-left
          x_spline_end,  0;          % 2: bottom-right = spline lower end
          x_spline_end,  Ly;         % 3: top-right    = spline upper end
          0,             Ly];        % 4: top-left

%% ===================== STEP III. SHAPE FUNCTIONS =======================
kWave = omega_eff / waveVelocity;   
nWaveFunctionsSet1 = 2*n1 + 2;
nWaveFunctionsSet2 = 2*n2 + 2;
nWaveFunctions     = nWaveFunctionsSet1 + nWaveFunctionsSet2;

% --- III.1 Wave numbers (two sets) --------------------------------------
% waves varying along x, propagating/evanescent in ±y-direction
kxy1(2*n1+2,2) = 0; %vector that stores the wave numbers
for n = 0:n1
    % ------- evaluate wave numbers and store them in a vector kxy1  ------
    kx1 = n*pi/Lx;
    ky1 = sqrt(kWave^2 - (n*pi/Lx)^2);
    kxy1(2*n+1,:) = [ kx1,  ky1];
    kxy1(2*n+2,:) = [ kx1, -ky1];
end

% waves varying along y, propagating/evanescent in ±x-direction
kxy2(2*n2+2,2) = 0; %vector that stores the wave numbers
for n = 0:n2
    % ------- evaluate wave numbers and store them in a vector kxy2  ------
    ky2 = n*pi/Ly;
    kx2 = sqrt(kWave^2 - (n*pi/Ly)^2);
    kxy2(2*n+1,:) = [ kx2,  ky2];
    kxy2(2*n+2,:) = [-kx2,  ky2];
end
kxy = [kxy1; kxy2];

%% =========== STEP IV. STIFFNESS MATRIX (boundary integrals) ============

% ---- Normals for straight edges (right edge uses pointwise normal) -----
tVecEdges = zeros(4,2);
nVecEdges = zeros(4,2);
for ii = 1:4
    if ii ~= 2      % edges 1,3,4: straight segments
        startPoint        = nodes(ii,:);
        endPoint          = nodes(mod(ii,4)+1,:);
        % Tangential vector along the straight edge
        tVecEdges(ii,:)   = endPoint - startPoint;
        cLtVec            = norm(tVecEdges(ii,:));
        % Outward normal (for a rectangular domain, orientation is standard)
        nVecEdges(ii,:)   = [ tVecEdges(ii,2), -tVecEdges(ii,1) ] / cLtVec;
    else
        % Edge 2 is the B-spline
        % its normal will be computed pointwise from the spline tangent
        nVecEdges(ii,:)   = [NaN, NaN];
    end
end

% ---- Gauss points & weights on each edge -------------------------------
xGP = zeros(4,nGP);
yGP = zeros(4,nGP);
wGP = zeros(4,nGP);
nVecCurve = [];           % size 2×nGP normals

for ii = 1:4
    if ii ~= 2
        % Straight edges
        startPoint = nodes(ii,:);
        endPoint   = nodes(mod(ii,4)+1,:);
        [x_,y_,w_] = getGPcoordinates(startPoint,endPoint,nGP);
        xGP(ii,:)  = x_';
        yGP(ii,:)  = y_';
        wGP(ii,:)  = w_';
    else
        % B-spline edge: GL integration, u ∈ [0,1]
        % with mapping to physical coordinates C(u) and Jacobian |C'(u)|
        [xhi, wStd] = getGauss1D(nGP);   % nodes/weights on [-1,1]

        % Map Gauss nodes/weights from [-1,1] to [0,1]:
        % u = (ξ+1)/2  →  w = (1/2) wStd   (Jacobian du/dξ = 1/2)
        uGL = 0.5*(xhi + 1);             
        wGL = 0.5*wStd;                  

        uGL = uGL.';                     % row vectors
        wGL = wGL.';

        nVecCurve  = zeros(2, nGP);
        for j = 1:nGP
            u = uGL(j);
            [C, dCdu] = evalBspline(u, B_control_points, B_degree, B_knot_vector);
            % [point, tangent] = evalBspline(...)
            xGP(ii,j) = C(1);
            yGP(ii,j) = C(2);

            % Tangent vector and arc-length factor for ds = |C'(u)| du
            % Line integral: ∫_Γ f ds = ∫_0^1 f(C(u)) |C'(u)| du
            % w_j,Γ = w_j,u · |C'(u_j)|
            t  = dCdu;
            Lt = norm(t);
            wGP(ii,j) = wGL(j) * Lt;   % physical line element weight

            % Outward normal
            nVecCurve(:,j) = -[ t(2); -t(1) ] / Lt;
        end
    end
end

K = zeros(nWaveFunctions, nWaveFunctions);

% ---------------- IV.1 Double loop over wave functions ------------------
tic;
for ii = 1:nWaveFunctions
    if ii <= 2*n1+2, set_ii = 1; else, set_ii = 2; end

    for jj = 1:nWaveFunctions
        if jj <= 2*n1+2, set_jj = 1; else, set_jj = 2; end

        % --------------- IV.2 Sum contributions from all edges ----------
        % Boundary integral over the entire closed boundary, decomposed
        % into 3 straight edges + 1 B-spline edge
        for kk = 1:4
            cBC  = BCs(kk);
            cWGP = wGP(kk,:);
            cXGP = xGP(kk,:);
            cYGP = yGP(kk,:);

            if kk == 2
                n1p = nVecCurve(1,:);
                n2p = nVecCurve(2,:);
            else
                n1p = nVecEdges(kk,1)*ones(1,nGP);
                n2p = nVecEdges(kk,2)*ones(1,nGP);
            end

            switch cBC
                case 'v'    % rigid wall: ∂p/∂n = 0  → v_n = 0
                    Psi_ii = evalShapeFunction(kxy(ii,1),kxy(ii,2),Lx,Ly,cXGP,cYGP,set_ii);
                    Psidot_jj_vec = evalShapeFunctionDerivative(kxy(jj,1),kxy(jj,2),Lx,Ly,cXGP,cYGP,set_jj);
                    Psidot_jj = Psidot_jj_vec(1,:).*n1p + Psidot_jj_vec(2,:).*n2p;

                    K(ii,jj)  = K(ii,jj) + (1i/(density*omega_eff)) * (cWGP.*Psi_ii)*Psidot_jj.';

                case 'p'    % pressure-release: p = 0
                    Psidot_ii_vec = evalShapeFunctionDerivative(kxy(ii,1),kxy(ii,2),Lx,Ly,cXGP,cYGP,set_ii);
                    Psidot_ii = Psidot_ii_vec(1,:).*n1p + Psidot_ii_vec(2,:).*n2p;

                    Psi_jj = evalShapeFunction(kxy(jj,1),kxy(jj,2),Lx,Ly,cXGP,cYGP,set_jj);

                    K(ii,jj) = K(ii,jj) ...
                               + (-1i/(density*omega_eff)) * (cWGP.*Psidot_ii) * Psi_jj.';
            end
        end
    end
end
toc;

%% ================== STEP V. LOAD VECTOR (boundary terms) ===============
f = zeros(nWaveFunctions,1);

for ii = 1:nWaveFunctions
    if ii <= 2*n1+2, set_ii = 1; else, set_ii = 2; end

    for kk = 1:4
        cBC  = BCs(kk);
        cWGP = wGP(kk,:);
        cXGP = xGP(kk,:);
        cYGP = yGP(kk,:);

        if kk == 2
            n1p = nVecCurve(1,:);
            n2p = nVecCurve(2,:);
        else
            n1p = nVecEdges(kk,1)*ones(1,nGP);
            n2p = nVecEdges(kk,2)*ones(1,nGP);
        end

        switch cBC
            case 'v'
                % Rigid wall
                Psi_ii = evalShapeFunction(kxy(ii,1),kxy(ii,2),Lx,Ly,cXGP,cYGP,set_ii);
                v_q_vec = evalLoadFunctionDerivative(q0,xQ,yQ,density,omega_eff,kWave,cXGP,cYGP);
                v_q = v_q_vec(1,:).*n1p + v_q_vec(2,:).*n2p;

                f(ii) = f(ii) - (cWGP.*Psi_ii)*v_q.';

            case 'p'
                % Pressure-release
                Psidot_ii_vec = evalShapeFunctionDerivative(kxy(ii,1),kxy(ii,2),Lx,Ly,cXGP,cYGP,set_ii);
                Psidot_ii = Psidot_ii_vec(1,:).*n1p + Psidot_ii_vec(2,:).*n2p;

                p_q = evalLoadFunction(q0,xQ,yQ,density,omega_eff,kWave,cXGP,cYGP);

                f(ii) = f(ii) + (1i/(density*omega_eff)) * (cWGP.*Psidot_ii) * p_q.';
        end

    end
end

%% ================== STEP VI. SOLVE ====================================
% Slight Tikhonov regularization to improve conditioning (as often used
% in WBM to deal with irregular frequencies, cf. Deckers 2010)
eps_reg = 1e-10;
w = (K + eps_reg*eye(size(K))) \ f;

%% ================== STEP VII. POST ====================================
% Resolution of the visualization grid
% Increasing N gives smoother and higher-resolution plots,
% but too large N makes the 3-D surfaces dense and harder to view.
% Adjust as needed
N = 120; 
[xSol, ySol] = meshgrid(linspace(0,Lx,N), linspace(0,Ly,N));
p_total = zeros(size(xSol));

% --- VII.1 Homogeneous part --------------------------------------------
for ii = 1:nWaveFunctions
    if ii <= 2*n1+2
        set_ii = 1;
    else
        set_ii = 2;
    end
    Psi_ii  = evalShapeFunction(kxy(ii,1),kxy(ii,2),Lx,Ly,xSol,ySol,set_ii);
    p_total = p_total + Psi_ii * w(ii);
end

% --- VII.2 Particular part ---------------------------------------------
p_part  = evalLoadFunction(q0,xQ,yQ,density,omega_eff,kWave,xSol,ySol);
p_total = p_total + p_part;

% === VII.3 Build B-spline curve (for geometry plot & masking) ===========
u_curve_plot = linspace(0,1,400);
B_curve = zeros(numel(u_curve_plot),2);
for i = 1:numel(u_curve_plot)
    [C,~] = evalBspline(u_curve_plot(i), B_control_points, B_degree, B_knot_vector);
    B_curve(i,:) = C;
end

% --- Geometry plot (to verify boundary vs. bounding box) ----------------
figure('Name','Domain Geometry');
hold on; axis equal; grid on;

% Three straight edges
plot([0, x_spline_end],[0,0],'b-','LineWidth',2);          
plot([0,0],[0,Ly],'b-','LineWidth',2);                     
plot([x_spline_end,0],[Ly,Ly],'b-','LineWidth',2);         

% Outer rectangular bounding box (original box used by WBM basis)
plot([x_spline_end,Lx],[0,0],'k--','LineWidth',1);
plot([x_spline_end,Lx],[Ly,Ly],'k--','LineWidth',1);
plot([Lx,Lx],[0,Ly],'k--','LineWidth',1);

% B-spline curve + control points + Gauss points + outward normals
plot(B_curve(:,1), B_curve(:,2), 'r-','LineWidth',2.5);    % spline
plot(B_control_points(:,1),B_control_points(:,2),'ro','MarkerFaceColor','y');
plot(xGP(2,:), yGP(2,:), 'go','MarkerSize',4,'MarkerFaceColor','g'); % GP
normal_scale = 0.07;
quiver(xGP(2,:), yGP(2,:), ...
       normal_scale*nVecCurve(1,:), normal_scale*nVecCurve(2,:), ...
       0, 'm', 'LineWidth', 1.2);

% Point source
plot(xQ,yQ,'c*','MarkerSize',10,'LineWidth',2);            

xlabel('x [m]'); ylabel('y [m]');
legend('cavity bottom','cavity left','cavity top', ...
       'box (bottom/right/top)','','', ...
       'B-spline','ctrl pts','Gauss pts','normals','source', ...
       'Location','best');
hold off;

% === VII.4 Use B-spline to build a mask (keep interior of cavity) =======

% Sort B-spline points by y to build x_bdry(y)
[y_sorted, idx] = sort(B_curve(:,2));
x_sorted = B_curve(idx,1);
y_min = min(y_sorted);  
y_max = max(y_sorted);

% For each grid point (xSol,ySol), interpolate the boundary x_bdry at the
% same horizontal line y.
y_query = ySol;
y_query(y_query < y_min) = y_min;
y_query(y_query > y_max) = y_max;
x_bdry_grid = interp1(y_sorted, x_sorted, y_query, 'linear');

% mask_outside = true  → points outside the cavity (will be set to NaN)
mask_outside = xSol > x_bdry_grid;

% === VII.5 Normalize and split Re / Im, apply mask ======================
% Use only interior points of the cavity to compute the maximum amplitude.
max_int = max(abs(p_total(~mask_outside)), [], 'all');
if isempty(max_int) || max_int == 0
    max_int = 1;
end

% Normalize such that max(|p|) ≈ 1 inside the cavity
p_norm = p_total / max_int;

pR = real(p_norm);
pI = imag(p_norm);

% Apply mask: everything outside the cavity becomes NaN
pR(mask_outside) = NaN;
pI(mask_outside) = NaN;

% === VII.6 3D Re / Im (with height, normalized, masked) =================
zscale = 1;   % Increase to accentuate height in 3D view

figure('Name','Re(p) normalized (3D, masked)');
surf(xSol, ySol, pR * zscale);
shading faceted;
axis tight;
xlabel('x [m]'); ylabel('y [m]'); zlabel('Re(p) (norm.)');
title(sprintf('Re(p) (normalized) at %.0f Hz',frequency));
colorbar; grid on; box on; view(3);

figure('Name','Im(p) normalized (3D, masked)');
surf(xSol, ySol, pI * zscale);
shading faceted;
axis tight;
xlabel('x [m]'); ylabel('y [m]'); zlabel('Im(p) (norm.)');
title(sprintf('Im(p) (normalized) at %.0f Hz',frequency));
colorbar; grid on; box on; view(3);

% === VII.7 |p| 2D contour (normalized, masked) ==========================
p_abs_norm = abs(p_norm);
p_abs_norm(mask_outside) = NaN;

figure('Name','|p| (normalized) contour, masked');
contourf(xSol, ySol, p_abs_norm, 30, 'LineStyle','none');
axis equal tight; colorbar;
xlabel('x [m]'); ylabel('y [m]');
title(sprintf('|p| (normalized) at %.0f Hz',frequency));

% === VII.8 |p| 2D contour (un-normalized, masked) =======================
p_abs_raw = abs(p_total);   % raw values (without division by max_int)
p_abs_raw(mask_outside) = NaN;

figure('Name','|p| (un-normalized) contour, masked');
contourf(xSol, ySol, p_abs_raw, 50, 'LineStyle','none');
axis equal tight;
colormap(jet);
colorbar;
xlabel('x [m]'); ylabel('y [m]');
title(sprintf('|p| (un-normalized) at %.0f Hz', frequency));

uiwait(gcf);

%% ======================== HELPER FUNCTIONS =============================

function Psi = evalShapeFunction(kx,ky,Lx,Ly,meshX,meshY,set)
% Shape functions for the two wave sets used in the WBM on a rectangle.

if set==1
    scalingFactor_y = double(imag(ky)>0);
    Psi = (cos(kx*meshX)).*exp(-1i*ky*(meshY - scalingFactor_y*Ly));
else
    scalingFactor_x = double(imag(kx)>0);
    Psi = exp(-1i*kx*(meshX - scalingFactor_x*Lx)).*(cos(ky*meshY));
end
end

function Psi_dot_vec = evalShapeFunctionDerivative(kx,ky,Lx,Ly,meshX,meshY,set)
% Gradient of the shape functions with respect to x and y.

if set==1
    scalingFactor_y = double(imag(ky)>0);
    Psi_dot_x = -kx * exp(-1i*ky*(meshY - scalingFactor_y*Ly)) .* sin(kx*meshX);
    Psi_dot_y = -1i*ky * exp(-1i*ky*(meshY - scalingFactor_y*Ly)) .* cos(kx*meshX);
else
    scalingFactor_x = double(imag(kx)>0);
    Psi_dot_x = -1i*kx * exp(-1i*kx*(meshX - scalingFactor_x*Lx)) .* cos(ky*meshY);
    Psi_dot_y = -ky    * exp(-1i*kx*(meshX - scalingFactor_x*Lx)) .* sin(ky*meshY);
end
Psi_dot_vec = [Psi_dot_x; Psi_dot_y];
end

function p_q = evalLoadFunction(p0,xQ,yQ,rho0,omega_eff,kWave,meshX,meshY)
% Particular solution for the sound pressure of a 2D point source
% (Hankel function of the 2nd kind, order 0)

r = hypot(meshX - xQ, meshY - yQ);
H_02 = besselj(0, kWave*r) - 1i*bessely(0, kWave*r);
p_q  = p0 * (rho0*omega_eff)/4 * H_02;
end

function v_q = evalLoadFunctionDerivative(p0,xQ,yQ,rho0,omega_eff,kWave,meshX,meshY)
% Gradient of the particular solution (for the normal velocity)

r = hypot(meshX - xQ, meshY - yQ);
r(r==0) = eps; % avoid division by zero
dHdx = -kWave*besselj(1,kWave*r).*((meshX-xQ)./r) + 1i*kWave*bessely(1,kWave*r).*((meshX-xQ)./r);
dHdy = -kWave*besselj(1,kWave*r).*((meshY-yQ)./r) + 1i*kWave*bessely(1,kWave*r).*((meshY-yQ)./r);
v_q  = +1i/(rho0*omega_eff)*p0*(rho0*omega_eff)/4 * [dHdx; dHdy];
end

function [xGP,yGP,wGP] = getGPcoordinates(startPoint,endPoint,nGP)
% getGPcoordinates: GL points on a line segment mapped
% from [-1,1] to the segment in global coordinates

% Obtain nodes/weights on [-1,1] (tabulated GL values)
[xhi, wGP_] = getGauss1D(nGP);  

% Length of the edge
lengthEdge = hypot(endPoint(1)-startPoint(1), endPoint(2)-startPoint(2));

% Map from [-1,1] to the actual line segment
xGP = ( endPoint(1) + startPoint(1) )/2 + xhi * ( endPoint(1)-startPoint(1) )/2;
yGP = ( endPoint(2) + startPoint(2) )/2 + xhi * ( endPoint(2)-startPoint(2) )/2;

jacobiDet = lengthEdge/2;
wGP = wGP_ * jacobiDet;

xGP = xGP.';
yGP = yGP.';
wGP = wGP.';
end

function [xhi, wGP_] = getGauss1D(nGP)
% getGauss1D: return Gauss–Legendre nodes xhi and weights wGP_ on [-1,1]
% for given nGP, using standard tabulated values

switch nGP
    case 1
        xhi  =   0 ;
        wGP_ =   2 ;
    case 2
        xhi  = [-0.577350269189626;0.577350269189626];
        wGP_ = [ 1.00000000000000;1.00000000000000];
    case 3
        xhi  = [-0.774596669241484;0;0.774596669241483];
        wGP_ = [ 0.555555555555556;0.888888888888889;0.555555555555556];
    case 4
        xhi  = [-0.861136311594053;-0.339981043584856;0.339981043584856;0.861136311594053];
        wGP_ = [ 0.347854845137454;0.652145154862546;0.652145154862546;0.347854845137454];
    case 5
        xhi  = [-0.906179845938664;-0.538469310105683;4.75059282543674e-17;0.538469310105683;0.906179845938664];
        wGP_ = [0.236926885056189;0.478628670499367;0.568888888888889;0.478628670499367;0.236926885056189];
    case 6
        xhi  = [-0.932469514203152;-0.661209386466265;-0.238619186083197;0.238619186083197;0.661209386466265;0.932469514203152];
        wGP_ = [0.171324492379171;0.360761573048139;0.467913934572691;0.467913934572692;0.360761573048139;0.171324492379170];
    case 7
        xhi  = [-0.949107912342758;-0.741531185599395;-0.405845151377397;-1.15122081766977e-16;0.405845151377397;0.741531185599395;0.949107912342759];
        wGP_ = [0.129484966168870;0.279705391489277;0.381830050505119;0.417959183673469;0.381830050505119;0.279705391489276;0.129484966168870];
    case 8
        xhi  = [-0.960289856497536;-0.796666477413627;-0.525532409916329;-0.183434642495650;0.183434642495650;0.525532409916329;0.796666477413627;0.960289856497536];
        wGP_ = [0.101228536290376;0.222381034453374;0.313706645877887;0.362683783378362;0.362683783378362;0.313706645877887;0.222381034453374;0.101228536290376];
    case 9
        xhi  = [-0.968160239507626;-0.836031107326636;-0.613371432700590;-0.324253423403809;9.31221438354534e-17;0.324253423403809;0.613371432700591;0.836031107326636;0.968160239507626];
        wGP_ = [0.0812743883615745;0.180648160694858;0.260610696402935;0.312347077040002;0.330239355001259;0.312347077040002;0.260610696402936;0.180648160694857;0.0812743883615745];
    case 10        
        xhi  = [-0.973906528517171;-0.865063366688985;-0.679409568299025;-0.433395394129247;-0.148874338981631;0.148874338981631;0.433395394129247;0.679409568299024;0.865063366688984;0.973906528517172];             
        wGP_ = [ 0.0666713443086883;0.149451349150581;0.219086362515981;0.269266719309996;0.295524224714752;0.295524224714753;0.269266719309996;0.219086362515982;0.149451349150581;0.0666713443086884];
    case 11
        xhi  = [-0.978228658146057;-0.887062599768096;-0.730152005574050;-0.519096129206812;-0.269543155952345;-1.56429560276319e-16;0.269543155952345;0.519096129206812;0.730152005574049;0.887062599768095;0.978228658146057];
        wGP_ = [0.0556685671161737;0.125580369464905;0.186290210927734;0.233193764591991;0.262804544510247;0.272925086777901;0.262804544510247;0.233193764591991;0.186290210927734;0.125580369464905;0.0556685671161738];  
    case 15
        xhi  = [-0.987992518020486;-0.937273392400706;-0.848206583410427;-0.724417731360170;-0.570972172608539;-0.394151347077563;-0.201194093997435;-1.04266455326068e-16;0.201194093997434;0.394151347077563;0.570972172608539;0.724417731360170;0.848206583410427;0.937273392400706;0.987992518020486];
        wGP_ = [0.0307532419961168;0.0703660474881088;0.107159220467172;0.139570677926155;0.166269205816994;0.186161000015562;0.198431485327112;0.202578241925561;0.198431485327111;0.186161000015562;0.166269205816994;0.139570677926154;0.107159220467172;0.0703660474881088;0.0307532419961167];
    case 50
        xhi  = [-0.998866404420071;-0.994031969432091;-0.985354084048006;-0.972864385106692;-0.956610955242808;-0.936656618944878;-0.913078556655792;-0.885967979523613;-0.855429769429946;-0.821582070859336;-0.784555832900399;-0.744494302226069;-0.701552468706822;-0.655896465685439;-0.607702927184950;-0.557158304514650;-0.504458144907464;-0.449806334974039;-0.393414311897565;-0.335500245419437;-0.276288193779532;-0.216007236876042;-0.154890589998146;-0.0931747015600862;-0.0310983383271890;0.0310983383271889;0.0931747015600862;0.154890589998146;0.216007236876042;0.276288193779532;0.335500245419437;0.393414311897565;0.449806334974039;0.504458144907464;0.557158304514650;0.607702927184951;0.655896465685440;0.701552468706822;0.744494302226068;0.784555832900399;0.821582070859336;0.855429769429946;0.885967979523613;0.913078556655792;0.936656618944878;0.956610955242808;0.972864385106692;0.985354084048006;0.994031969432091;0.998866404420071];
        wGP_ = [0.00290862255315519;0.00675979919574557;0.0105905483836510;0.0143808227614850;0.0181155607134894;0.0217802431701250;0.0253606735700125;0.0288429935805350;0.0322137282235786;0.0354598356151461;0.0385687566125877;0.0415284630901475;0.0443275043388029;0.0469550513039493;0.0494009384494663;0.0516557030695811;0.0537106218889967;0.0555577448062126;0.0571899256477284;0.0586008498132225;0.0597850587042660;0.0607379708417703;0.0614558995903161;0.0619360674206835;0.0621766166553467;0.0621766166553467;0.0619360674206835;0.0614558995903164;0.0607379708417699;0.0597850587042655;0.0586008498132226;0.0571899256477290;0.0555577448062123;0.0537106218889961;0.0516557030695811;0.0494009384494661;0.0469550513039493;0.0443275043388037;0.0415284630901482;0.0385687566125871;0.0354598356151455;0.0322137282235779;0.0288429935805351;0.0253606735700131;0.0217802431701251;0.0181155607134890;0.0143808227614853;0.0105905483836513;0.00675979919574539;0.00290862255315491];
    case 100
        xhi  = [-0.999713726773442;-0.998491950639596;-0.996295134733125;-0.993124937037444;-0.988984395242992;-0.983877540706057;-0.977809358486918;-0.970785775763706;-0.962813654255816;-0.953900782925492;-0.944055870136256;-0.933288535043080;-0.921609298145334;-0.909029570982530;-0.895561644970727;-0.881218679385018;-0.866014688497165;-0.849964527879591;-0.833083879888401;-0.815389238339176;-0.796897892390314;-0.777627909649496;-0.757598118519707;-0.736828089802021;-0.715338117573056;-0.693149199355802;-0.670283015603141;-0.646761908514129;-0.622608860203708;-0.597847470247179;-0.572501932621381;-0.546597012065094;-0.520158019881763;-0.493210789208191;-0.465781649773358;-0.437897402172032;-0.409585291678302;-0.380872981624630;-0.351788526372422;-0.322360343900529;-0.292617188038472;-0.262588120371504;-0.232302481844974;-0.201789864095736;-0.171080080538603;-0.140203137236114;-0.109189203580061;-0.0780685828134366;-0.0468716824215917;-0.0156289844215430;0.0156289844215431;0.0468716824215915;0.0780685828134368;0.109189203580061;0.140203137236114;0.171080080538603;0.201789864095736;0.232302481844974;0.262588120371504;0.292617188038472;0.322360343900529;0.351788526372422;0.380872981624630;0.409585291678302;0.437897402172032;0.465781649773358;0.493210789208191;0.520158019881763;0.546597012065094;0.572501932621381;0.597847470247179;0.622608860203708;0.646761908514129;0.670283015603141;0.693149199355802;0.715338117573056;0.736828089802021;0.757598118519707;0.777627909649495;0.796897892390315;0.815389238339176;0.833083879888401;0.849964527879591;0.866014688497165;0.881218679385019;0.895561644970727;0.909029570982530;0.921609298145334;0.933288535043079;0.944055870136256;0.953900782925492;0.962813654255816;0.970785775763707;0.977809358486918;0.983877540706057;0.988984395242992;0.993124937037444;0.996295134733125;0.998491950639596;0.999713726773442];
        wGP_ = [0.000734634490505451;0.00170939265351806;0.00268392537155344;0.00365596120132667;0.00462445006342261;0.00558842800386551;0.00654694845084481;0.00749907325546452;0.00844387146966869;0.00938041965369451;0.0103078025748689;0.0112251140231861;0.0121314576629796;0.0130259478929717;0.0139077107037186;0.0147758845274413;0.0156296210775463;0.0164680861761454;0.0172904605683238;0.0180959407221280;0.0188837396133749;0.0196530874944355;0.0204032326462093;0.0211334421125271;0.0218430024162471;0.0225312202563363;0.0231974231852546;0.0238409602659687;0.0244612027079572;0.0250575444815795;0.0256294029102078;0.0261762192395460;0.0266974591835711;0.0271926134465773;0.0276611982207923;0.0281027556591010;0.0285168543223952;0.0289030896011252;0.0292610841106383;0.0295904880599129;0.0298909795933328;0.0301622651051695;0.0304040795264549;0.0306161865839805;0.0307983790311522;0.0309504788504908;0.0310723374275667;0.0311638356962101;0.0312248842548492;0.0312554234538630;0.0312554234538630;0.0312248842548491;0.0311638356962100;0.0310723374275670;0.0309504788504906;0.0307983790311522;0.0306161865839805;0.0304040795264546;0.0301622651051693;0.0298909795933329;0.0295904880599124;0.0292610841106385;0.0289030896011252;0.0285168543223953;0.0281027556591013;0.0276611982207921;0.0271926134465769;0.0266974591835708;0.0261762192395457;0.0256294029102076;0.0250575444815797;0.0244612027079573;0.0238409602659685;0.0231974231852545;0.0225312202563364;0.0218430024162480;0.0211334421125277;0.0204032326462097;0.0196530874944352;0.0188837396133742;0.0180959407221280;0.0172904605683231;0.0164680861761453;0.0156296210775462;0.0147758845274410;0.0139077107037189;0.0130259478929718;0.0121314576629800;0.0112251140231860;0.0103078025748691;0.00938041965369427;0.00844387146966850;0.00749907325546483;0.00654694845084544;0.00558842800386549;0.00462445006342231;0.00365596120132615;0.00268392537155359;0.00170939265351769;0.000734634490505757];
    otherwise
        error('nGP=%d not implemented in getGauss1D.', nGP);
end
end

function [point, tangent] = evalBspline(u, control_points, degree, knot_vector)
% Evaluate a B-spline curve and its derivative at parameter u
% using Cox–de Boor recurrence (basis + derivative).

n = size(control_points,1) - 1;  % last index
point   = [0,0];
tangent = [0,0];
for i = 0:n
    N  = evalBsplineBasis(i, degree, u, knot_vector);
    dN = evalBsplineBasisDerivative(i, degree, u, knot_vector);
    point   = point   + N  * control_points(i+1,:);
    tangent = tangent + dN * control_points(i+1,:);
end
end

function N = evalBsplineBasis(i, p, u, U)
% B-spline basis function N_{i,p}(u) (Cox–de Boor)

if p==0
    if (u>=U(i+1) && u<U(i+2)) || (abs(u-U(end))<1e-12 && abs(U(i+2)-U(end))<1e-12)
        N = 1;
    else
        N = 0;
    end
else
    left = 0; right = 0;
    denL = U(i+p+1) - U(i+1);
    denR = U(i+p+2) - U(i+2);
    if denL > 1e-14  % if denL != 0
        left = (u - U(i+1))/denL * evalBsplineBasis(i,  p-1, u, U);
    end
    if denR > 1e-14 % if denR != 0
        right = (U(i+p+2)-u)/denR * evalBsplineBasis(i+1,p-1, u, U);
    end
    N = left + right;
end
end

function dN = evalBsplineBasisDerivative(i, p, u, U)
% Derivative of the B-spline basis function N_{i,p}(u).

if p==0
    dN = 0;
else
    denL = U(i+p+1) - U(i+1);
    denR = U(i+p+2) - U(i+2);
    left = 0; right = 0;
    if denL > 1e-14
        left  =  p/denL * evalBsplineBasis(i,  p-1, u, U);
    end
    if denR > 1e-14
        right = -p/denR * evalBsplineBasis(i+1,p-1, u, U);
    end
    dN = left + right;
end
end

