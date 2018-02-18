function [ output ] = xor_long(input_bits, xor_bit)

if xor_bit==0
    output = input_bits;
else 
    l = size(input_bits,2);
    output = mod(input_bits+ones(1,l),2);
end


end

