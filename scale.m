function [out_data_i, out_data_q] = scale(fi, fd_lo, t, in_data)

out_data_i = sign(cos(2*pi*(fi+fd_lo)*t)).* in_data;
out_data_q = sign(sin(2*pi*(fi+fd_lo)*t)).* in_data;

end