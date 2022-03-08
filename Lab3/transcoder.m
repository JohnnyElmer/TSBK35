function [bits, bpp, dist, psnr]=transcoder(im1, blocksize,qy,qc, transformType, huffmanOrJpeg, sampleFactor )

% This is a very simple transform coder and decoder. Copy it to your directory
% and edit it to suit your needs.
% You probably want to supply the image and coding parameters as
% arguments to the function instead of having them hardcoded.


% Read an image
im=double(imread(im1))/255;

% What blocksize do we want?
%Sent in from main.

% Quantization steps for luminance and chrominance
% Sent in from main.

% Change colourspace 
imy=rgb2ycbcr(im);

bits=0;

% Somewhere to put the decoded image
imyr=zeros(size(im));

% First we code the luminance component
% Here comes the coding part

if strcmp(transformType, 'dct')
    tmp = bdct(imy(:,:,1), blocksize); % DCT

elseif strcmp(transformType, 'dwht')
    tmp = bdwht(imy(:,:,1), blocksize); % DWHT   

end

tmp = bquant(tmp, qy);             % Simple quantization

if strcmp(huffmanOrJpeg, 'huffman')
    p = ihist(tmp(:));                 % Only one huffman code
    bits = bits + huffman(p);          % Add the contribution from
                                       % each component
elseif strcmp(huffmanOrJpeg, 'jpeg')
    bits = sum(jpgrate(tmp, blocksize));

elseif strcmp(huffmanOrJpeg, 'sepHuffman')
    for k= 1:size(tmp,1)
        p = ihist(tmp(k,:));
        bits = bits + huffman(p);
    end
end

% Here comes the decoding part
tmp = brec(tmp, qy);               % Reconstruction

if strcmp(transformType, 'dct')
    imyr(:,:,1) = ibdct(tmp, blocksize, [512 768]);  % Inverse DCT

elseif strcmp(transformType, 'dwht')
    imyr(:,:,1) = ibdwht(tmp, blocksize, [512 768]);  % Inverse DWHT

end

% Next we code the chrominance components
for c=2:3                          % Loop over the two chrominance components
  % Here comes the coding part

  tmp = imy(:,:,c);

  % If you're using chrominance subsampling, it should be done
  % here, before the transform.
  if (sampleFactor ~= 0)
  tmp = imresize(tmp, 1 / sampleFactor);
  chromaSize = size(tmp);
  end

  if strcmp(transformType, 'dct')
    tmp = bdct(tmp, blocksize); % DCT

  elseif strcmp(transformType, 'dwht')
    tmp = bdwht(tmp, blocksize); %DWHT
  end
  
  tmp = bquant(tmp, qc);           % Simple quantization
  
  if strcmp(huffmanOrJpeg, 'huffman')
    p = ihist(tmp(:));                 % Only one huffman code
    bits = bits + huffman(p);          % Add the contribution from
                                       % each component
  elseif strcmp(huffmanOrJpeg, 'jpeg')
        bits = bits + sum(jpgrate(tmp, blocksize));

  elseif strcmp(huffmanOrJpeg, 'ind huffman')
    for k=1:size(tmp, 1)
        p = ihist(tmp(k, :));
        bits = bits + huffman(p);
    end
    
  end
			
  % Here comes the decoding part
  tmp = brec(tmp, qc);            % Reconstruction
  
  if strcmp(transformType, 'dct')
      if (sampleFactor == 0)
      tmp = ibdct(tmp, blocksize, [512 768]);  % Inverse DCT
        
      elseif (sampleFactor ~= 0)
        tmp = ibdct(tmp, blocksize, chromaSize);  % Inverse DCT
      end

  elseif strcmp(transformType, 'dwht')
       if (sampleFactor == 0)
            tmp = ibdwht(tmp, blocksize, [512 768]); % Inverse DWHT
       elseif (sampleFactor ~= 0)
            tmp = ibdwht(tmp, blocksize, chromaSize); % Inverse DWHT
       end
  end

  % If you're using chrominance subsampling, this is where the
  % signal should be upsampled, after the inverse transform.
   if (sampleFactor ~= 0)
        tmp = imresize(tmp, sampleFactor);
   end
  imyr(:,:,c) = tmp;
  
end

% Display total number of bits and bits per pixel
bits;
bpp = bits/(size(im,1)*size(im,2));

% Revert to RGB colour space again.
imr=ycbcr2rgb(imyr);

% Measure distortion and PSNR
dist = mean((im(:)-imr(:)).^2);
psnr = 10*log10(1/dist);

% Display the original image
%figure, imshow(im)
%title('Original image')

%Display the coded and decoded image
%figure, imshow(imr);
%title(sprintf('Decoded image, %5.2f bits/pixel, PSNR %5.2f dB', bpp, psnr))

