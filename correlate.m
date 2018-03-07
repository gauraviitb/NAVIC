function [acq_mat_ele] = correlate (in_phase, Q_comp, incoming)
i = 0;
j = 0;
%i = corrcoef(in_phase,incoming);
%j = corrcoef(Q_comp,incoming);
L = size(in_phase);%
for V = 1: L

i = i + circshift(in_phase, [V, 0]) * incoming' ;
j = j + circshift(Q_comp, [V, 0]) * incoming' ;

end
acq_mat_ele = i*i + j*j ;
if(i<0)
acq_mat_ele = -1 * acq_mat_ele;
end
end