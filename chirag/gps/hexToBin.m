function [ bin_array] = hexToBin( str )

m = length(str);
bin_array = zeros(1,4*m);

for i =0:m-1
    a = [0 0 0 0];
    switch str(i+1)
        case '0'
            a = [0 0 0 0];
        case '1'
            a = [0 0 0 1];
        case '2'
            a = [0 0 1 0];
        case '3'
            a = [0 0 1 1];
        case '4'
            a = [0 1 0 0];
        case '5'
            a = [0 1 0 1];
        case '6'
            a = [0 1 1 0];
        case '7'
            a = [0 1 1 1];
        case '8'
            a = [1 0 0 0];
        case '9'
            a = [1 0 0 1];
        case 'A'
            a = [1 0 1 0];
        case 'B'
            a = [1 0 1 1];
        case 'C'
            a = [1 1 0 0];
        case 'D'
            a = [1 1 0 1];
        case 'E'          
            a = [1 1 1 0];
        case 'F'
            a = [1 1 1 1];
        otherwise
            disp('Error')
            break
    end
     bin_array(4*i+1:4*(i+1))= a;       
end

end

