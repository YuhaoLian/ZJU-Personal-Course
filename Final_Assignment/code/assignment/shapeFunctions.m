function [N, dN] = shapeFunctions(L1, L2, L3)
    N = [L1*(2*L1-1);
         L2*(2*L2-1);
         L3*(2*L3-1);
         4*L1*L2;
         4*L2*L3;
         4*L3*L1];
    
    % 形函数导数 (对L1和L2)
    dN = zeros(6,2);
    dN(1,:) = [4*L1-1, 0];             % dN1/dL1, dN1/dL2
    dN(2,:) = [0, 4*L2-1];             % dN2/dL1, dN2/dL2
    dN(3,:) = [0, 0];                  % dN3/dL1, dN3/dL2 (将用链式法则处理)
    dN(4,:) = [4*L2, 4*L1];            % dN4/dL1, dN4/dL2
    dN(5,:) = [-4*L2, 4*(1-2*L2-L1)];  % dN5/dL1, dN5/dL2
    dN(6,:) = [4*(1-2*L1-L2), -4*L1];  % dN6/dL1, dN6/dL2
    
    % 处理dN3 (使用链式法则)
    dN3_dL3 = 4*L3 - 1;
    dN(3,:) = [-dN3_dL3, -dN3_dL3];  % dL3/dL1 = -1, dL3/dL2 = -1
end