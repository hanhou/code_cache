function cm = colormap_BlueWhiteRed(n, gamma, offset)

if nargin<1; n = 100; end
if nargin<2; gamma = 0.6; end
if nargin<3; offset = 0; end

n_pos = round((1 - offset) * n);
n_neg = round((1 + offset) * n);

% cm = ([n*ones(1,n), n:-1:0 ; ...
%       0:n, n-1:-1:0; ...
%       0:n, ones(1,n)*n]' / n).^gamma;
%   

c_pos = ([n_pos*ones(1,n_pos); 0:n_pos-1; 0:n_pos-1]' / n_pos).^gamma;
c_neg = ([n_neg-1:-1:0; n_neg-1:-1:0; ones(1,n_neg)*n_neg]' / n_neg).^gamma;

cm = [c_pos ; 1,1,1; c_neg];
cm = cm(end:-1:1,:);  

colormap(gca, cm);
