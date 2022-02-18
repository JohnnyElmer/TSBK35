function x=imdct(y)
%imdct - Inverse MDCT of transformed signal
%------------------------------------------------------------------------------
%SYNOPSIS	X = imdct(Y)
%		  Perform blockwise IMDCT on a transformed matrix Y.
%                 Each column is Y is treated as a transform block and
%                 the window size W is double the column length of Y.
%                 The function assumes that data has been padded with
%                 zeros at each end before being transformed, so these
%                 zeros are removed after the inverse transform.
%                 The return value X is a column vector.
%
%                 The window used is sin(pi/W*(k+1/2)).
%
%
%NOTE           The transform is not normalized.
%
%
%SEE ALSO	mdct
%
%------------------------------------------------------------------------------
%Harald Nautsch                        (C) 2007 Image Coding Group. LiU, SWEDEN

%RCSID          $Id: imdct.m,v 1.2 2007/02/16 15:22:53 harna Exp $

if (nargin == 0)
  error('No input arguments.')
end



ws=2*size(y, 1);
  
w=sin(pi*((0:ws-1)+0.5)/ws)';

T=zeros(ws,ws/2);

for m=1:ws/2
  T(:, m) = cos(pi/(2*ws)*((0:2:2*ws-2)+1+ws/2)*(2*m-1))';
end

x = T*y;

x = x.*repmat(w, 1, size(y,2));

x = 4/ws * x;

x = reshape(x, ws/2, []);

x = x(:, 2:end-1);

x = x(:, 1:2:end)+x(:, 2:2:end);

x=x(:);

