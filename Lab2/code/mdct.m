function y=mdct(x, ws)
%mdct - MDCT of signal
%------------------------------------------------------------------------------
%SYNOPSIS	Y = mdct(X, W)
%		  Perform blockwise MDCT on a signal vector X, using
%                 a window size W. The resulting transformed data Y
%                 is a matrix with W/2 rows, one column for each block.
%                 The data is padded with W/2 zeros at each end before
%                 being transformed.
%                 The window used is sin(pi/W*(k+1/2)).
%
%
%NOTE           The transform is not normalized.
%
%
%SEE ALSO	imdct
%
%------------------------------------------------------------------------------
%Harald Nautsch                        (C) 2007 Image Coding Group. LiU, SWEDEN

%RCSID          $Id: mdct.m,v 1.2 2007/02/16 15:22:59 harna Exp $

if (nargin == 0)
  error('No input arguments.')
end

if (nargin == 1)
  error('No window size given.')
end


x=x(:);
  
y=zeros(ws, 2*length(x)/ws+1);

y(:,1) = [zeros(ws/2,1); x(1:ws/2)];
y(:,end) = [x(end-ws/2+1:end); zeros(ws/2,1)];

y(:,2:2:end) = reshape(x,ws,[]);
y(:,3:2:end-1) = reshape(x(1+ws/2:end-ws/2),ws,[]);

w=sin(pi*((0:ws-1)+0.5)/ws)';

y=y.*repmat(w, 1, size(y,2));

T=zeros(ws/2,ws);

for m=1:ws/2
  T(m, :) = cos(pi/(2*ws)*((0:2:2*ws-2)+1+ws/2)*(2*m-1));
end

y=T*y;
