% Created: Nov, 07, 2025 11:23:01 by Wataru Fukuda
clc; clear; close all;

%%% configuration %%%
nodes = [ 0, 0;
          1, 2;
          2, 0;
          3, 2];  % control points
p = 2;            % polynomial degree
division = 10;    % subdivision per element

%%% data preparation %%%
nn = size(nodes,1); % number of nodes
if p+1 > nn
  error("Number of nodes should be greater than (polynomial degree + 1)")
end
ne = nn - p; % number of elements
kv = makeKnotVector(p,nn,ne);
unique_kv = unique(kv);
xi_all = [];
for e = 1:ne
  xi_start = unique_kv(e);
  xi_end = unique_kv(e+1);
  xi_all = [xi_all, linspace(xi_start, xi_end, division+1)];
end
xi_all = unique(xi_all);

%%% evaluate B-spline curve %%%
curve = zeros(length(xi_all),2);
for k = 1:length(xi_all)
  xi = xi_all(k);
  N = zeros(1,nn);
  for i = 1:nn
    N(i) = bsplineBasis(i,p,kv,xi);
  end
  curve(k,:) = N*nodes;
end

%%% visualization %%%
figure; hold on; grid on; axis equal;
plot(nodes(:,1), nodes(:,2), 'ro--', 'LineWidth', 1.5);
plot(curve(:,1), curve(:,2), 'b-', 'LineWidth', 2);
legend('Control Points', 'B-spline Curve');
xlabel('x'); ylabel('y'); title('B-spline Curve Visualization');

%%% keep figures %%%
uiwait(gcf);

function kv = makeKnotVector(p,nn,ne);
  kv = [zeros(1,p+1), 1:(ne-1), ne*ones(1,p+1)];
  % kv = kv/max(kv); % normalization, not necessary
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


