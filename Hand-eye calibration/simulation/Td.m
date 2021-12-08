function [T] = Td(x)
%平移变换方程


T  = [1         0         0        x(1)
      0         1         0        x(2)
      0         0         1        x(3)
      0         0         0        1];
  

end

