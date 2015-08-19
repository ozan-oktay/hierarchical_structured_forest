function [ output_img ] = imRescale( input_img, min_new, max_new, save_type )
%IMRESCALE Summary of this function goes here
%   Detailed explanation goes here

if (nargin<4)
  save_type = class(input_img);
end

max_val  = max(input_img(:));
min_val  = min(input_img(:));

slope      = (max_new-min_new) / (max_val-min_val);
intersec   = min_new - slope * min_val;
input_img  = input_img * slope + intersec;
output_img = cast(input_img,save_type);

end

