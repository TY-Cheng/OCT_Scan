function [ x ] = ladexp_huberTV( z,B,par )
%PROXIMAL Takes a noisy OCT image outputs an estimate for the latent image
%
%   Inputs:
%
%       z: noisy OCT image (dimension: r*c)
%       B: PSF (convolution kernel), if no blur, B=1
%
%       par:
%       par.gamma: (1/2)weight on F in prox_{\gamma F} (default: 2)
%       par.tau: (1/2)weight on G in prox_{prox_{\tau G}} (default: 0.5)
%       par.lambda: weight on TV regularization
%       par.theta: extrapolation stepsize (default: 0.3)
%       par.c1: distribution mean coefficient (default: 0.959)
%       par.c2: distribution variance coefficient (default: 0.0804)
%       par.maxIter: maximum number of iterations (default: 10)
%
%   Output:
%
%       x: estimated latent image (r*c)
%
%-------------------------------------------------------------

%-------------------------------------------------------------
%   Parameter selection
%-------------------------------------------------------------
% default parameters
par0.gamma = 0.5;
par0.tau = 0.5;
par0.lambda = 0.5;
par0.theta = 0.3;
par0.c1 = 0.959;
par0.c2 = 0.0804;
par0.avrg = 1;
par0.maxIter = 30;

if exist('par','var')
    if ~isstruct(par)
        error('dualrecons: parameter structure is not a matlab struct');
    end
    % merge default with given values
    params = mergeParam(par0, par);
else
    params = par0;
end

lambda = params.lambda;
c1 = params.c1;
c2 = params.c2;
avrg = params.avrg;
maxIter = params.maxIter;

if ~isfloat(z)
    z=im2double(z);
end

%-------------------------------------------------------------
%   Initialization
%-------------------------------------------------------------
z(isnan(z)) = 1e-4;
xold = log(conv2(z,fspecial('gaussian',[3 3], 0.4),'same')); %%% averaging is important for avoiding generating fake high values of x
% xold = log(z);
xbar = xold;
[m,n,c] = size(xbar);
yold = zeros(m,2*n);
pold = zeros(m,2*n);

alpha = 2;
delta = 0.0018;
gamma = 0.98;
%%% prior operators

dxf=[1 -1];
dyf=[1;-1];

% K1 = @(u) [conv2([u,u(:,end)],dxf,'valid'), conv2([u;u(end,:)],dyf,'valid')];
% K1_T = @(hy) rot90(conv2([rot90(hy(:,1:n),2), rot90(hy(:,1),2)],dxf,'valid'),2) + ...
%             rot90(conv2([rot90(hy(:,n+1:2*n),2); rot90(hy(1,n+1:2*n),2)],dyf,'valid'),2);        

yold(:,1:2*n) = K1(xold);

%-------------------------------------------------------------
%   Linearized ADMM on exponential form with huber TV
%-------------------------------------------------------------

for iter = 1:maxIter
    %%% update x
    df =  -1/(2*c2)*z.*exp(-xold) + c1/(2*c2)*sqrt(z).*exp(-xold/2) + (1-0.5*avrg);
    x = xold - delta*( df - alpha*K1_T((yold-K1(xold))) - K1_T(pold)  );
    
    K1x = K1(x);
    %%% update y
    thresh_huber = 0.02; % based on the range of image values
    range = max(max(K1x)) - min(min(K1x));
    pd = (K1x - pold/alpha)/range;
    pdnorm = sqrt( pd(:,1:n).^2 + pd(:,n+1:2*n).^2 );
%     AA = zeros(size(pdnorm));
%     AA(pdnorm-lambda/alpha > pdnorm/(1+lambda/(alpha*thresh_huber))) = 1;
%     figure; imshow(AA);
    y = range* repmat(max(pdnorm-lambda/alpha, pdnorm/(1+lambda/(alpha*thresh_huber))),1,2).*(pd./repmat(max(pdnorm,1e-5),1,2));
%     v = K1(x)-pold/alpha;
%     vnorm = sqrt( v(:,1:n).^2 + v(:,n+1:2*n).^2 );
%     y = repmat(max(vnorm - lambda/alpha, 0)./max(vnorm,1e-6), 1,2) .* v;
    
    %%% update p
    p = pold + alpha*(y - K1x);
    
    
    %%% correction step
    dx = (1/delta)*(xold-x) - df +  -1/(2*c2)*z.*exp(-x) + c1/(2*c2)*sqrt(z).*exp(-x/2)+(1-0.5*avrg);
    dy = alpha*(yold-y) - alpha*K1(xold-x);
    dp = (1/alpha)*(pold-p);
    r = norm(dx,'fro')^2 + norm(dy,'fro')^2 + norm(dp,'fro')^2;
    phi = reshape((pold-p),1,[])*reshape((K1(xold)-yold),[],1) + ...
        reshape((xold-x),1,[])*reshape(dx,[],1) + ...
        reshape((yold-y),1,[])*reshape(dy,[],1);
    beta = gamma*phi/r;
    
    x = x-beta*dx;
    y = y-beta*dy;
    p = p-beta*dp;
    
    %%% stopping
    if sum(reshape((x-xold).^2,[],1))/sqrt(numel(x)) <0.0005
       break
    else
        yold = y;
        xold = x;
        pold = p;

        %%% compute energy
%         E = lambda*sum(sum(abs(K1(x)))) + sum(sum(1/(2*c2)*( sqrt(z./x)-c1 ).^2 + 0.5*log(x)))
%         E = lambda*sum(sum(abs(K1(x)))) + sum(sum(1/(2*c2)*( sqrt(z).*exp(-x/2)-c1 ).^2 + (1-0.5*avrg)*x))
        %%% display images
        

    end
    
end

x = exp(x);

end

function [x] = K1(u)
    x = [u(:,[2:end end])-u, u([2:end end],:,:)-u];
end

function [u] = K1_T(x)
    n = size(x,2)/2;
    u = ( -x(:,1:n) + x(:,[1 1:n-1]) ) + ( -x(:,(n+1):2*n) + x([1 1:end-1],(n+1):2*n) );
end
