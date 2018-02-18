function [out_data] = shift(in_data)
out_data = circshift(in_data', 1) ;
out_data = out_data';
end;