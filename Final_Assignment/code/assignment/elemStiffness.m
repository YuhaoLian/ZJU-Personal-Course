function [Ke, Fe] = elemStiffness(vertices)
    % 高斯积分点（面积坐标）- 7点积分
    gauss_points = [
        1/3, 1/3, 1/3;
        0.0597158717, 0.4701420641, 0.4701420641;
        0.4701420641, 0.0597158717, 0.4701420641;
        0.4701420641, 0.4701420641, 0.0597158717;
        0.7974269853, 0.1012865073, 0.1012865073;
        0.1012865073, 0.7974269853, 0.1012865073;
        0.1012865073, 0.1012865073, 0.7974269853
    ];
    
    weights = [
        0.225;
        0.1323941527;
        0.1323941527;
        0.1323941527;
        0.1259391805;
        0.1259391805;
        0.1259391805
    ];
    
    Ke = zeros(6,6);
    Fe = zeros(6,1);
    
    for i = 1:7
        L1 = gauss_points(i, 1);
        L2 = gauss_points(i, 2);
        L3 = gauss_points(i, 3);
        
        % 计算形函数和导数
        [N, dN] = shapeFunctions(L1, L2, L3);
        
        % 计算雅可比矩阵 (修正为2x2)
        dxdL = [dN(:,1)'*vertices(:,1), dN(:,1)'*vertices(:,2);
                dN(:,2)'*vertices(:,1), dN(:,2)'*vertices(:,2)];
        
        J = dxdL;  % 2x2雅可比矩阵
        detJ = abs(det(J));
        
        % ========== 关键修改3：避免奇异矩阵 ==========
        if detJ < 1e-14
            error('雅可比矩阵奇异，单元退化');
        end
        
        invJ = inv(J);
        
        % 转换导数到物理坐标系
        dNdx = zeros(6,2);
        for j = 1:6
            dNdl = [dN(j,1); dN(j,2)];
            dNdx(j,:) = (invJ * dNdl)';
        end
        
        % 计算积分点物理坐标
        x_phys = N' * vertices(:,1);
        y_phys = N' * vertices(:,2);
        
        % 计算Ke项
        for m = 1:6
            for n = 1:6
                Ke(m,n) = Ke(m,n) + (dNdx(m,:)*dNdx(n,:)' + 4*N(m)*N(n)) * detJ * weights(i);
            end
            % 计算Fe项
            Fe(m) = Fe(m) + f([x_phys, y_phys]) * N(m) * detJ * weights(i);
        end
    end
end