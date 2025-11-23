% Created: Nov, 10, 2025 16:10:44 by Wataru Fukuda
%% ########################################################################
%  #-------------------- STEP I - PREPARING PROGRAM ----------------------#
%  ########################################################################

clear classes; clear all; clc


%% ########################################################################
%  #-------------------- STEP II - DEFINE INPUT DATA ---------------------#
%  ########################################################################

% config = 'init';
config = 'Bspline';

switch config
  case 'init'
    nodes = { [ 0       0 ;
                20      0 ] ;
              [ 20      0 ;
                20      10 ] ;
              [ 20      10 ;
                0       10 ] ;
              [ 0       10 ;
                0       0 ]
            } ;
    Lx = 20;
    Ly = 10;
        
    n1 = 30;
    n2 = 30;

    nGP = 50;

    q0 = 1;
    xQ = 2;
    yQ = 4;

    BCs = ['v';
           'p';
           'v';
           'p' ] ;

    density      = 1.25;
    waveVelocity = 343;
    frequency    = 50;               % frequency in Hz
    omega        = frequency*2*pi;   % circular frequency

  case 'Bspline'
    nodes = { [ 0   0 ;
                1   0 ] ;
              [ 1   0 ;
                1.16  0.16 ;
                1.34  0.50 ;
                1.16  0.84 ;
                1   1 ];
              [ 1   1 ;
                0   1 ] ;
              [ 0   1 ;
                0   0 ]
            } ;
    p = 3;

    Lx = 1.25;  % todo
    Ly = 1;

    n1 = 5;  % todo
    n2 = 5;

    nGP = 15;  % todo

    q0 = 1;
    xQ = 0.25;
    yQ = 0.25;

    BCs = ['v';
           'v';
           'v';
           'v' ] ;

    density      = 1.225;
    waveVelocity = 339;
    frequency    = 800;               % frequency in Hz
    omega        = frequency*2*pi;   % circular frequency

end


%% ########################################################################
%  #-------------------- STEP III - SHAPE FUNCTIONS ----------------------#
%  ########################################################################

kWave = omega/waveVelocity;
nWaveFunctionsSet1 = 2*n1+2;
nWaveFunctionsSet2 = 2*n2+2;
nWaveFunctions = nWaveFunctionsSet1+nWaveFunctionsSet2;

%  ************************************
%  ***  III.1 Evaluate Wave Numbers ***
%  ************************************

kxy1(nWaveFunctionsSet1,2) = 0;         %vector that stores the wave numbers

for n=0:1:n1  % eq10
  kx1 = n*pi/Lx;
  ky1 = sqrt(kWave^2-kx1^2);
  kxy1(2*n+1,1) =  kx1;
  kxy1(2*n+1,2) =  ky1;
  kxy1(2*n+2,1) =  kx1;
  kxy1(2*n+2,2) = -ky1;
end

kxy2(nWaveFunctionsSet2,2) = 0;         %vector that stores the wave numbers

for n=0:1:n2  % eq11
  ky2 = n*pi/Ly;
  kx2 = sqrt(kWave^2-ky2^2);
  kxy2(2*n+1,1) = kx2;
  kxy2(2*n+1,2) = ky2;
  kxy2(2*n+2,1) = -kx2;
  kxy2(2*n+2,2) = ky2;
end

kxy = [kxy1;kxy2];


%% ########################################################################
%  #--------------- STEP IV - CALCULATE STIFFNESS MATRIX -----------------#
%  ########################################################################

xGP = zeros(4,nGP);
yGP = zeros(4,nGP);
wGP = zeros(4,nGP);
nVecEdges = cell(4,1);

for ii = 1:4
  [xhi, w] = getGP(nGP);
  if size(nodes{ii}) == 2
    startPoint = nodes{ii}(1,:);
    endPoint = nodes{ii}(2,:);
    [xGP_,yGP_,wGP_] = calcGPcoordinates(startPoint,endPoint,xhi,w);
    xGP(ii,:) = xGP_';
    yGP(ii,:) = yGP_';
    wGP(ii,:) = wGP_';
    nVecEdges{ii} = calcNormalVectors(startPoint, endPoint);
  else
    [xGP_,yGP_,wGP_] = calcGPcoordinatesBspline(nodes{ii},p,xhi,w);
    xGP(ii,:) = xGP_';
    yGP(ii,:) = yGP_';
    wGP(ii,:) = wGP_';
    nVecEdges{ii} = calcNormalVectorsBspline(nodes{ii}, p, xhi, w);
  end
end

K = zeros(nWaveFunctions,nWaveFunctions);

%  ***********************************************************
%  ***  IV.1 Loop to calculate entries of stiffness matrix ***
%  ***********************************************************
%        _    _  
%    /\   |    |   dim1: Matrix Dimension equals to number of 
%     dim1  | A ii(v)|     Element DOFs
%    \/   |_    _| 
%     
%         <- dim1 ->  
%
% size(A ii): |2*(n1+n2)+4| x |2*(n1+n2)+4|


%  *************************************
%  ***  Loop over Shape Functions ii ***
%  ************************************* 

tic

for ii = 1: nWaveFunctions
  
  % if ii <= 2*n1+2
  if ii <= nWaveFunctionsSet1
    set_ii = 1;
  else 
    set_ii = 2;
  end
  
  %  *************************************
  %  ***  Loop over Shape Functions jj ***
  %  *************************************  
  
  for jj = 1: nWaveFunctions

    % if jj <= 2*n1+2
    if jj <= nWaveFunctionsSet1
      set_jj = 1;
    else 
      set_jj = 2;
    end
      
    %  *********************************
    %  ***  IV.2 Loop over all edges ***
    %  *********************************
    %           ____4____
    %           /      | 
    %  		       1 /       |3
    %         /      |
    %        /_______2_____|		
    %		       	
    %  calculate stiffness contribution for each edge separately
    %        
    %        /           |
    %    k=1  /    + k=2   + k=3 |   + ...
    %      /             |
    %       /   _______2_____    |

    for kk = 1:4
      
      cBC  = BCs(kk);
      cWGP = wGP(kk,:);
      cXGP = xGP(kk,:);
      cYGP = yGP(kk,:);
      nVec = nVecEdges{kk};
      
      switch cBC

        case 'v'
          
          % evaluate Shape Function 
          % at the Gauß Point coordinates
          Psi_ii = evalShapeFunction(kxy(ii,1),kxy(ii,2),...
                 Lx,Ly,...
                 cXGP,cYGP,set_ii);
               
          % evaluate Shape Function Derivatives 
          % as vector quantity
          % at the Gauß Point coordinates          
          Psidot_jj_vec = evalShapeFunctionDerivative...
                (kxy(jj,1),kxy(jj,2),...
                 Lx,Ly,...
                 cXGP,cYGP,set_jj);
          % calculate normal derivative of sound velocity
          % in outwards direction perpendicular to the edge
          Psidot_jj =  Psidot_jj_vec(1,:) * nVec(1) + ...
                 Psidot_jj_vec(2,:) * nVec(2);
                         
          % carrying out the numerical integration
          % here done in form of a vector scalar product
          
          K(ii,jj)  = (1i/(density*omega)* ...
                  cWGP.*Psi_ii)*Psidot_jj.' + K(ii,jj);  % todo: where does this equation come from?
              
        case 'p'
          
          % evaluate Shape Function Derivatives 
          % as vector quantity
          % at the Gauß Point coordinates          
          Psidot_ii_vec = evalShapeFunctionDerivative...
                          (kxy(ii,1),kxy(ii,2),...
                           Lx,Ly,...
                           cXGP,cYGP,set_ii);  
          % calculate normal derivative of sound velocity
          % in outwards direction perpendicular to the edge
          Psidot_ii =  Psidot_ii_vec(1,:) * nVec(1) + ...
                 Psidot_ii_vec(2,:) * nVec(2);
               
          % evaluate Shape Function 
          % at the Gauß Point coordinates                         
          Psi_jj = evalShapeFunction(kxy(jj,1),kxy(jj,2),...
                           Lx,Ly,...
                           cXGP,cYGP,set_jj);
                         
          % carrying out the numerical integration
          % here done in form of a vector scalar product
          
          K(ii,jj) = (-1i/(density*omega)* ...
                  cWGP.*Psidot_ii)*Psi_jj.' + K(ii,jj);  % todo: where does this equation come from?
              
        case 'z'
          K(ii,jj) = 0;
      
          
      end

    end

  end
  
end

toc

% Optional - Matrix Calculation instead of Loop
% Uncoment block to use


%% ########################################################################
%  #------------------ STEP V - CALCULATE LOAD VECTOR --------------------#
%  ########################################################################

f = zeros(nWaveFunctions,1);

%  *****************************************************
%  ***  V.1 Loop to calculate entries of load vector ***
%  *****************************************************

for ii = 1:nWaveFunctions
  
  % if ii <= 2*n1+2
  if ii <= nWaveFunctionsSet1
    set_ii = 1;
  else 
    set_ii = 2;
  end  
    
%  ********************************
%  ***  V.2 Loop over all edges ***
%  ********************************
  for kk = 1:4
    
    cBC  = BCs(kk);
    cWGP = wGP(kk,:);
    cXGP = xGP(kk,:);
    cYGP = yGP(kk,:);
    nVec = nVecEdges{kk};
    
    switch cBC

      case 'v'
        
        % evaluate Shape Function 
        % at the Gauß Point coordinates
        Psi_ii = evalShapeFunction(kxy(ii,1),kxy(ii,2),...
               Lx,Ly,...
               cXGP,cYGP,set_ii);
             
        % evaluate load function derivative at Gauß Points     
        v_q_vec  = evalLoadFunctionDerivative(q0,xQ,yQ,density,omega,kWave,cXGP,cYGP);   
        % calculate normal derivative of sound velocity
        % in outwards direction perpendicular to the edge
        v_q =  v_q_vec(1,:) * nVec(1) + ...
             v_q_vec(2,:) * nVec(2); 
             
        f(ii) = - (cWGP.*Psi_ii)*v_q.' + f(ii);  % todo: where does this equation come from?
      case 'p'
        
        % evaluate Shape Function Derivatives 
        % as vector quantity
        % at the Gauß Point coordinates          
        Psidot_ii_vec = evalShapeFunctionDerivative...
                        (kxy(ii,1),kxy(ii,2),...
                         Lx,Ly,...
                         cXGP,cYGP,set_ii);  
        % calculate normal derivative of sound velocity
        % in outwards direction perpendicular to the edge
        Psidot_ii =  Psidot_ii_vec(1,:) * nVec(1) + ...
               Psidot_ii_vec(2,:) * nVec(2);  
        % evaluate load function at Gauß Points     
        p_q   = evalLoadFunction(q0,xQ,yQ,density,omega,kWave,cXGP,cYGP);
        
        % carrying out the numerical integration
        % here done in form of a vector scalar product          
        f(ii) = (1i/(density*omega)* ...
                 cWGP.*Psidot_ii)*p_q.' + f(ii);  % todo: where does this equation come from?
      case 'z'
      f(ii) = 0;
        
    end 
    
  end
end

%% ########################################################################
%  #---------------- STEP VI - SOLVE SYSTEM OF EQUATIONS -----------------#
%  ########################################################################

w = linsolve(K,f);

%% ########################################################################
%  #---------------- STEP VII - RESUME AND PLOT RESULTS -----------------#
%  ########################################################################

% define mesh for plot
switch config
  case "init"
    [xSol,ySol] = meshgrid(linspace(0,Lx,50),linspace(0,Ly,50));
  case "Bspline"
    etaLine = linspace(-1,1,50);
    [Xi,Eta] = meshgrid(etaLine, etaLine);
    xL = zeros(size(etaLine));
    yL = (etaLine+1)/2;
    [xR,yR,wGP_] = calcGPcoordinatesBspline(nodes{2},p,etaLine.',w);
    xL = repmat(xL(:), 1, length(etaLine));
    yL = repmat(yL(:), 1, length(etaLine));
    xR = repmat(xR(:), 1, length(etaLine));
    yR = repmat(yR(:), 1, length(etaLine));
    xSol = (1 - Xi) / 2 .* xL + (Xi + 1) / 2 .* xR;
    ySol = (1 - Xi) / 2 .* yL + (Xi + 1) / 2 .* yR;
end

p_total   = zeros(size(xSol));

%  **********************************************
%  ***  VII.1 Resume Homogenous solution part ***
%  **********************************************

% Homogenous solution
for ii = 1:nWaveFunctions
  
  % if ii <= 2*n1+2
  if ii <= nWaveFunctionsSet1
    set_ii = 1;
  else 
    set_ii = 2;
  end 
  
  Psi_ii  = evalShapeFunction(kxy(ii,1),kxy(ii,2),...
        Lx,Ly,...
        xSol,ySol,set_ii); 
      
  p_total = Psi_ii * w(ii) + p_total;

end

%  *******************************************
%  ***  VII.2 Add Particular solution part ***
%  *******************************************

% Particular solution
p_part  = evalLoadFunction(q0,xQ,yQ,density,omega,kWave,xSol,ySol);
p_total = p_part + p_total;

%  ********************************
%  ***  VII.3 Visualize results ***
%  ********************************

% plot real part of total solution
% figure('Name','Real Part of Total Solution','NumberTitle','off');
% surf(xSol,ySol,real(p_total))
%
% % plot imaginary part of total solution
% figure('Name','Imaginary Part of Total Solution','NumberTitle','off');
% surf(xSol,ySol,imag(p_total))
%
% % plot a single shape function
% figure('Name','Single Shape Function','NumberTitle','off');
% ax1 = subplot(2,1,1);
% surf(xSol,ySol,real(evalShapeFunction(kxy(3,1),kxy(3,2),...
%         Lx,Ly,...
%         xSol,ySol,1)))
% ax1 = subplot(2,1,2);
% surf(xSol,ySol,imag(evalShapeFunction(kxy(3,1),kxy(3,2),...
%         Lx,Ly,...
%         xSol,ySol,1)))      

figure('Name','Total Solution Visualization','NumberTitle','off','Position',[100 100 1200 800]);

subplot(2,2,1);
surf(xSol, ySol, real(p_total)); title('Real Part of Total Solution');
xlabel('x'); ylabel('y'); shading interp; colorbar;

subplot(2,2,2);
surf(xSol, ySol, imag(p_total)); title('Imaginary Part of Total Solution');
xlabel('x'); ylabel('y'); shading interp; colorbar;

subplot(2,2,3);
surf(xSol, ySol, real(evalShapeFunction(kxy(3,1),kxy(3,2), Lx, Ly, xSol, ySol, 1)));
title('Real Part of Shape Function');
xlabel('x'); ylabel('y'); shading interp; colorbar;

subplot(2,2,4);
surf(xSol, ySol, imag(evalShapeFunction(kxy(3,1),kxy(3,2), Lx, Ly, xSol, ySol, 1)));
title('Imaginary Part of Shape Function');
xlabel('x'); ylabel('y'); shading interp; colorbar;

sgtitle('Total Solution and Single Shape Function');

%%% keep figures %%%
uiwait(gcf);

% #########################################################################
%             
%       FUNCTION: evalShapeFunction
%
% *************************************************************************
% 
% Input: kx, ky    -> wave numbers
%    Lx, Ly    -> maximal dimensions of the problem's circumscribing
%             rectangular
%    meshX, meshY  -> meshgrid of coordinates of GP along the edges,
%             x coordinate and belonging y coordinate of a
%             point are stored at the same position
%
% *************************************************************************
% 
% Description:  Wiki
%
% #########################################################################

function Psi = evalShapeFunction(kx,ky,Lx,Ly,meshX,meshY,set)

  if set == 1
    
    % ------- evaluate wave numbers and store them in a vector kxy1  ------
    
    scalingFactor_y = 0;
      
    if imag(ky) > 0
      scalingFactor_y = 1;
    end
    
    % --------- evaluate wave functions at all points of the mesh grid ----
    % --------- defined by matrices x and y              ----
    
    Psi = (cos(kx*meshX)).*exp(-1i*ky*(meshY-scalingFactor_y*Ly));  % todo: what is this scaling factor?

  elseif set == 2 % shapeFunction or DOF belongs to second set
    
    % ------- evaluate wave numbers and store them in a vector kxy2  ------
    
    scalingFactor_x = 0;
    if imag(kx) > 0
      scalingFactor_x = 1;
    end
    
    % --------- evaluate wave functions at all points of the mesh grid ----
    % --------- defined by matrices x and y              ----
    
    Psi = exp(-1i*kx*(meshX-scalingFactor_x*Lx)).*(cos(ky*meshY));
     
  end

end

% #########################################################################
%             
%       FUNCTION: evalShapeFunctionDerivatives
%
% *************************************************************************
% 
% Input: kx, ky    -> wave numbers
%    Lx, Ly    -> maximal dimensions of the problem's circumscribing
%             rectangular
%    meshX, meshY  -> meshgrid of coordinates of GP along the edges,
%             x coordinate and belonging y coordinate of a
%             point are stored at the same position
%
% *************************************************************************
% 
% Description:  Wiki
%
% #########################################################################

function Psi_dot_vec = evalShapeFunctionDerivative(kx,ky,Lx,Ly,meshX,meshY,set)

  if set == 1 
    
    scalingFactor_y = 0;
    
    if imag(ky)>0
      scalingFactor_y = 1;
    end
    
    % --------- evaluate wave functions at all points of the mesh grid ----
    % --------- defined by matrices x and y              ----
    
    Psi_dot_x = -kx * exp(-ky*(meshY-scalingFactor_y*Ly)*1i).* ...
           sin(kx*meshX);
    
    % --------- evaluate wave functions at all points of the mesh grid ----
    % --------- defined by matrices x and y              ----
    
    Psi_dot_y = -ky * exp(-ky*(meshY-scalingFactor_y*Ly)*1i).* ...
           cos(kx*(meshX))*1i;
    
  elseif set == 2
    
    scalingFactor_x = 0;
    if imag(kx)>0
      scalingFactor_x = 1;
    end
        
    % --------- evaluate wave functions at all points of the mesh grid ----
    % --------- defined by matrices x and y              ----
    
    Psi_dot_x = -kx * exp(-kx*(meshX-scalingFactor_x*Lx)*1i).* ...
           cos(ky*meshY)*1i;
      
    % --------- evaluate wave functions at all points of the mesh grid ----
    % --------- defined by matrices x and y              ----
    
    Psi_dot_y = -ky * exp(-kx*(meshX-scalingFactor_x*Lx)*1i).* ...
           sin(ky*meshY);

  end

Psi_dot_vec = [Psi_dot_x;Psi_dot_y];

end

% #########################################################################
%             
%       FUNCTION: evalLoadFunction
%
% *************************************************************************
% 
% Input: 
%    meshX, meshY  -> meshgrid of coordinates of GP along the edges,
%             x coordinate and belonging y coordinate of a
%             point are stored at the same position
%
% *************************************************************************
% 
% Description:  Wiki
%
% #########################################################################

function p_q = evalLoadFunction(p0,xQ,yQ,rho0,omega,kWave,meshX,meshY)

% -- set up Hankel function of 2nd kind and evaluate alongside all edge  --
  
  H_02 =  besselj(0, kWave*((meshX - xQ).^2 + (meshY - yQ).^2).^(1/2)) ...
      - 1i * bessely(0, kWave*((meshX - xQ).^2 + (meshY - yQ).^2).^(1/2));
    
% --------- particular solution for sound pressure and sound velocity -----
       
  p_q = p0 * (rho0*omega)/4 * H_02;  
  
end

% #########################################################################
%             
%       FUNCTION: evalLoadFunctionDerivative
%
% *************************************************************************
% 
% Input: 
%    meshX, meshY  -> meshgrid of coordinates of GP along the edges,
%             x coordinate and belonging y coordinate of a
%             point are stored at the same position
%
% *************************************************************************
% 
% Description:  Wiki
%
% #########################################################################

function v_q = evalLoadFunctionDerivative(p0,xQ,yQ,rho0,omega,kWave,meshX,meshY)

% -- set up derivative of Hankel function of 2nd kind and evaluate alongside all edges --

  dH_02_dx = - (kWave * besselj( 1, kWave*((meshX - xQ).^2 + ...
                     (meshY - yQ).^2).^(1/2)) .* ...
         (2*meshX - 2*xQ) ) ./ ...
         (2*((meshX - xQ).^2 + (meshY - yQ).^2).^(1/2)) ...
         + (kWave * bessely( 1, kWave*((meshX - xQ).^2 + ...
         (meshY - yQ).^2).^(1/2)) .* ...
         (2*meshX - 2*xQ)*1i) ./ ...
         (2*((meshX - xQ).^2 + (meshY - yQ).^2).^(1/2));

  dH_02_dy = - (kWave * besselj( 1, kWave *((meshX - xQ).^2 + ...
                      (meshY - yQ).^2).^(1/2)) .* ...
                (2*meshY - 2*yQ)) ./ ...
         (2*((meshX - xQ).^2 + (meshY - yQ).^2).^(1/2)) ...
         + (kWave * bessely(1, kWave *((meshX - xQ).^2 + ...
                      (meshY - yQ).^2).^(1/2)) .* ...
                (2*meshY - 2*yQ)*1i) ./ ...
         (2*((meshX - xQ).^2 + (meshY - yQ).^2).^(1/2));
       
% --------- particular solution for sound pressure and sound velocity -----

  v_q = +1i/(rho0*omega)*p0*(rho0*omega)/4 * [dH_02_dx ; dH_02_dy]; 
                         
end


function [xhi, wGP_] = getGP(nGP)
% xhi = -1..1
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
    case 12
    case 13
    case 14
    case 15
      xhi  = [-0.987992518020486;-0.937273392400706;-0.848206583410427;-0.724417731360170;-0.570972172608539;-0.394151347077563;-0.201194093997435;-1.04266455326068e-16;0.201194093997434;0.394151347077563;0.570972172608539;0.724417731360170;0.848206583410427;0.937273392400706;0.987992518020486];
      wGP_ = [0.0307532419961168;0.0703660474881088;0.107159220467172;0.139570677926155;0.166269205816994;0.186161000015562;0.198431485327112;0.202578241925561;0.198431485327111;0.186161000015562;0.166269205816994;0.139570677926154;0.107159220467172;0.0703660474881088;0.0307532419961167];
    case 50
      xhi  = [-0.998866404420071;-0.994031969432091;-0.985354084048006;-0.972864385106692;-0.956610955242808;-0.936656618944878;-0.913078556655792;-0.885967979523613;-0.855429769429946;-0.821582070859336;-0.784555832900399;-0.744494302226069;-0.701552468706822;-0.655896465685439;-0.607702927184950;-0.557158304514650;-0.504458144907464;-0.449806334974039;-0.393414311897565;-0.335500245419437;-0.276288193779532;-0.216007236876042;-0.154890589998146;-0.0931747015600862;-0.0310983383271890;0.0310983383271889;0.0931747015600862;0.154890589998146;0.216007236876042;0.276288193779532;0.335500245419437;0.393414311897565;0.449806334974039;0.504458144907464;0.557158304514650;0.607702927184951;0.655896465685440;0.701552468706822;0.744494302226068;0.784555832900399;0.821582070859336;0.855429769429946;0.885967979523613;0.913078556655792;0.936656618944878;0.956610955242808;0.972864385106692;0.985354084048006;0.994031969432091;0.998866404420071];
      wGP_ = [0.00290862255315519;0.00675979919574557;0.0105905483836510;0.0143808227614850;0.0181155607134894;0.0217802431701250;0.0253606735700125;0.0288429935805350;0.0322137282235786;0.0354598356151461;0.0385687566125877;0.0415284630901475;0.0443275043388029;0.0469550513039493;0.0494009384494663;0.0516557030695811;0.0537106218889967;0.0555577448062126;0.0571899256477284;0.0586008498132225;0.0597850587042660;0.0607379708417703;0.0614558995903161;0.0619360674206835;0.0621766166553467;0.0621766166553467;0.0619360674206835;0.0614558995903164;0.0607379708417699;0.0597850587042655;0.0586008498132226;0.0571899256477290;0.0555577448062123;0.0537106218889961;0.0516557030695811;0.0494009384494661;0.0469550513039493;0.0443275043388037;0.0415284630901482;0.0385687566125871;0.0354598356151455;0.0322137282235779;0.0288429935805351;0.0253606735700131;0.0217802431701251;0.0181155607134890;0.0143808227614853;0.0105905483836513;0.00675979919574539;0.00290862255315491];
    case 100
      xhi  = [-0.999713726773442;-0.998491950639596;-0.996295134733125;-0.993124937037444;-0.988984395242992;-0.983877540706057;-0.977809358486918;-0.970785775763706;-0.962813654255816;-0.953900782925492;-0.944055870136256;-0.933288535043080;-0.921609298145334;-0.909029570982530;-0.895561644970727;-0.881218679385018;-0.866014688497165;-0.849964527879591;-0.833083879888401;-0.815389238339176;-0.796897892390314;-0.777627909649496;-0.757598118519707;-0.736828089802021;-0.715338117573056;-0.693149199355802;-0.670283015603141;-0.646761908514129;-0.622608860203708;-0.597847470247179;-0.572501932621381;-0.546597012065094;-0.520158019881763;-0.493210789208191;-0.465781649773358;-0.437897402172032;-0.409585291678302;-0.380872981624630;-0.351788526372422;-0.322360343900529;-0.292617188038472;-0.262588120371504;-0.232302481844974;-0.201789864095736;-0.171080080538603;-0.140203137236114;-0.109189203580061;-0.0780685828134366;-0.0468716824215917;-0.0156289844215430;0.0156289844215431;0.0468716824215915;0.0780685828134368;0.109189203580061;0.140203137236114;0.171080080538603;0.201789864095736;0.232302481844974;0.262588120371504;0.292617188038472;0.322360343900529;0.351788526372422;0.380872981624630;0.409585291678302;0.437897402172032;0.465781649773358;0.493210789208191;0.520158019881763;0.546597012065094;0.572501932621381;0.597847470247179;0.622608860203708;0.646761908514129;0.670283015603141;0.693149199355802;0.715338117573056;0.736828089802021;0.757598118519707;0.777627909649495;0.796897892390315;0.815389238339176;0.833083879888401;0.849964527879591;0.866014688497165;0.881218679385019;0.895561644970727;0.909029570982530;0.921609298145334;0.933288535043079;0.944055870136256;0.953900782925492;0.962813654255816;0.970785775763707;0.977809358486918;0.983877540706057;0.988984395242992;0.993124937037444;0.996295134733125;0.998491950639596;0.999713726773442];
      wGP_ = [0.000734634490505451;0.00170939265351806;0.00268392537155344;0.00365596120132667;0.00462445006342261;0.00558842800386551;0.00654694845084481;0.00749907325546452;0.00844387146966869;0.00938041965369451;0.0103078025748689;0.0112251140231861;0.0121314576629796;0.0130259478929717;0.0139077107037186;0.0147758845274413;0.0156296210775463;0.0164680861761454;0.0172904605683238;0.0180959407221280;0.0188837396133749;0.0196530874944355;0.0204032326462093;0.0211334421125271;0.0218430024162471;0.0225312202563363;0.0231974231852546;0.0238409602659687;0.0244612027079572;0.0250575444815795;0.0256294029102078;0.0261762192395460;0.0266974591835711;0.0271926134465773;0.0276611982207923;0.0281027556591010;0.0285168543223952;0.0289030896011252;0.0292610841106383;0.0295904880599129;0.0298909795933328;0.0301622651051695;0.0304040795264549;0.0306161865839805;0.0307983790311522;0.0309504788504908;0.0310723374275667;0.0311638356962101;0.0312248842548492;0.0312554234538630;0.0312554234538630;0.0312248842548491;0.0311638356962100;0.0310723374275670;0.0309504788504906;0.0307983790311522;0.0306161865839805;0.0304040795264546;0.0301622651051693;0.0298909795933329;0.0295904880599124;0.0292610841106385;0.0289030896011252;0.0285168543223953;0.0281027556591013;0.0276611982207921;0.0271926134465769;0.0266974591835708;0.0261762192395457;0.0256294029102076;0.0250575444815797;0.0244612027079573;0.0238409602659685;0.0231974231852545;0.0225312202563364;0.0218430024162480;0.0211334421125277;0.0204032326462097;0.0196530874944352;0.0188837396133742;0.0180959407221280;0.0172904605683231;0.0164680861761453;0.0156296210775462;0.0147758845274410;0.0139077107037189;0.0130259478929718;0.0121314576629800;0.0112251140231860;0.0103078025748691;0.00938041965369427;0.00844387146966850;0.00749907325546483;0.00654694845084544;0.00558842800386549;0.00462445006342231;0.00365596120132615;0.00268392537155359;0.00170939265351769;0.000734634490505757];
      
  end
end

function [xGP_,yGP_,wGP_] = calcGPcoordinates(startPoint,endPoint,xhi,w)
  % transform from local convective coordinates to global coordinate system
  length = abs( sqrt( ( endPoint(1) - startPoint(1) )^2 + ( endPoint(2) - startPoint(2) )^2 ) );
  xGP_ = ( endPoint(1) + startPoint(1) ) / 2 + xhi * ( endPoint(1) - startPoint(1) ) / 2;
  yGP_ = ( endPoint(2) + startPoint(2) ) / 2 + xhi * ( endPoint(2) - startPoint(2) ) / 2;  
  jacobiDet = length/2;
  wGP_ = w * jacobiDet;
end

function [xGP_,yGP_,wGP_] = calcGPcoordinatesBspline(cps,p,xhi,w)
  nn = size(cps,1);
  ne = nn - p;
  kv = makeKnotVector(p,nn,ne);
  xhi = ( max(kv) + min(kv) ) / 2 + xhi * ( max(kv) - min(kv) ) / 2;
  coord = zeros(size(xhi,1),2);
  for xhi_index = 1:size(xhi,1)
    N = zeros(1,nn);
    for i = 1:nn
      N(i) = bsplineBasis(i,p,kv,xhi(xhi_index));
    coord(xhi_index,:) = N * cps;
    end
  end
  xGP_ = coord(:,1);
  yGP_ = coord(:,2);
  jacobiDet = 1;  % todo
  wGP_ = w * jacobiDet;
end

function nVec = calcNormalVectors(startPoint, endPoint)
  tVecEdges = zeros(1,2);
  tVecEdges = [endPoint(1)-startPoint(1) endPoint(2)-startPoint(2)];  % y-coordinate
  cLtVec      = abs(sqrt(tVecEdges(1)^2+tVecEdges(2)^2));
  nVec = [tVecEdges(2)/cLtVec,-tVecEdges(1)/cLtVec];
end

function nVecBspline = calcNormalVectorsBspline(cps, p, xhi, w)
  nn = size(cps,1);
  ne = nn - p; 
  kv = makeKnotVector(p,nn,ne);
  xhi = ( max(kv) + min(kv) ) / 2 + xhi * ( max(kv) - min(kv) ) / 2;
  nVecBspline = zeros(size(xhi,1),2);
  for j = 1:size(xhi,1)
    % tVec = [0,0];
    N = zeros(1,nn);
    for i = 1:nn
      N(i) = bsplineBasisDerivative(i,p,kv,xhi(j));
    end
    tVec = N * cps;
    cLtVec = abs(sqrt(tVec(1)^2+tVec(2)^2));
    nVec = [tVec(2)/cLtVec,-tVec(1)/cLtVec];
    nVecBspline(j,:) = nVec;
  end
end

function dN = bsplineBasisDerivative(i,p,kv,xi)
  if p==0
    dN = 0;
    return
  end
  a1 = kv(i+p) - kv(i);
  a2 = kv(i+p+1) - kv(i+1);
  left = 0; 
  right = 0;
  if a1 ~= 0
    left = p/a1 * bsplineBasis(i,p-1,kv,xi);
  end
  if a2 ~= 0
    right = p/a2 * bsplineBasis(i+1,p-1,kv,xi);
  end
  dN = left - right;
end

function N = bsplineBasis(i,p,kv,xi)
  if p == 0
  if kv(i) <= xi && xi < kv(i+1)
    N = 1;
  else
    N = 0;
  end
  if xi == kv(end) && xi == kv(i+1)
    N = 1;
  end
  else
  a1 = kv(i+p) - kv(i);
  a2 = kv(i+p+1) - kv(i+1);
  left = 0;
  right = 0;
  if a1 ~= 0
    left = (xi - kv(i)) / a1 * bsplineBasis(i,p-1,kv,xi);
  end
  if a2 ~= 0
    right = (kv(i+p+1) - xi) / a2 * bsplineBasis(i+1,p-1,kv,xi);
  end
  N = left + right;
  end
end

function kv = makeKnotVector(p,nn,ne);
  kv = [zeros(1,p+1), 1:(ne-1), ne*ones(1,p+1)];
  % kv = kv/max(kv); % normalization, not necessary
end
