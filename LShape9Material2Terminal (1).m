clc
clear
close all
tic;  % 开始计时
MP = struct( ...
    "Emax",   [320, 340, 360, 370, 380, 390, 400, 420, 510], ...  % GPa
    "density",[11.3,11.8,12.4,12.7,13.0,13.3,13.6,14.2,16.9], ...
    "cost",   [508.5,542.8,620.0,660.4,702.0,731.5,761.6,823.6,1284.4], ...
    "Emin",1e-3, ...
    "v",0.3, ...
    "P",1.5 ...
);

GP = struct("ngp",4,...% number of Gaussian points
            "noph",9,...% number of material phases
            "els",1,...% size of unit element
            "lx",60,"ly",60,"cutout_size",36);% geometry of structure
OP = struct("NOI",2500,...% number of iterations
            "Tho",2.5,...% 增加过滤强度，使结构更平滑
            "delta_T",0.01,...% 减小时间步长，提高稳定性
            "VC0",[0.15,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.15],...
            "Mu",2*ones(1,9), ...  % 减小惩罚系数
            "nRelax",10*ones(1,9), ...  % 减少松弛迭代次数
            "Gamma",0.28*ones(1,9),...  % 减小初始Gamma值
            "dGamma",0.005*ones(1,9), ...  % 减小Gamma增量
            "maxGamma",5*ones(1,9), ...  % 减小最大Gamma值
            "Lambda",0*ones(1,9), ...
            "PlotResult",1);
[ndof,noe,NodeCoord,EDofMat,ElNodes,NodeRep,Phi,IC,nelx,nely] = MeshGeneration(GP.lx,GP.ly,GP.els,GP.noph,GP.cutout_size);
[F,Dis,FixedNodes,ukdis,kdis] = BoundaryConditionsImplementation(ndof,NodeCoord);
iK = reshape(kron(EDofMat,ones(8,1))',8*8*noe,1);
jK = reshape(kron(EDofMat,ones(1,8))',8*8*noe,1);
[Ke,Vole,Bmat,Emat] = PlaneElementStiffnessCalculation(GP.ngp,GP.els,MP.v);
[T1e,T2e] = RD_T1T2Generator(GP.ngp,GP.els);
% T1 and T2 assembly
iT = reshape(kron(ElNodes,ones(4,1))',4*4*noe,1);
jT = reshape(kron(ElNodes,ones(1,4))',4*4*noe,1);
st1 = reshape(T1e(:)*ones(1,noe),4*4*noe,1);
st2 = reshape(T2e(:)*ones(1,noe),4*4*noe,1);
T1 = sparse(iT,jT,st1);
T2 = sparse(iT,jT,st2);
opc = 1;
[r,s] = meshgrid(linspace(-1,1,20));% [r,s] = Grid on each element for Phi interpolation
ns = numel(s);% ns is the number of points in each element for Phi interpolation
tmpPhi = zeros(numel(s),noe,GP.noph);
vlfe = zeros(GP.noph,noe);% vlfe = The percentage of materials from each phase in the element
VC = zeros(GP.noph,OP.NOI);
Compliance = zeros(OP.NOI,1)
vlfe_phi = zeros(GP.noph,noe,GP.noph);% derivative of volume fraction wrt phi
Lambda = zeros(GP.noph,1);% Lambda in augmented Lagrngian
while opc<=OP.NOI
    for i=1:GP.noph
        tmpPhi(:,:,i) = 0.25*((1-r(:)).*(1-s(:))*Phi(ElNodes(:,1),:,i)'+...
                              (1+r(:)).*(1-s(:))*Phi(ElNodes(:,2),:,i)'+...
                              (1-r(:)).*(1+s(:))*Phi(ElNodes(:,3),:,i)'+...% Phi =(结点数,1,noph)，tmpPhi = zeros(numel(s),noe,GP.noph);
                              (1+r(:)).*(1+s(:))*Phi(ElNodes(:,4),:,i)');% tmpPhi = Interpolated Phi values on element for volume fraction calculation
    end
    hphie = tmpPhi>0;% hphie = heavisided phi in elements
    for i=1:GP.noph
% Calculation of vlfe
        if i<GP.noph
            cond = sum(hphie(:,:,1:i),3)+(~hphie(:,:,i+1));
            vlfe(i,:) = sum(cond==(i+1))/ns;
            for j=1:i
                if (i==1 && j==1)
                    vlfe_phi(1,:,1) = sum(1-hphie(:,:,2))/ns;
                else
                    ind = [setdiff(1:i,j)];
                    cond = sum(hphie(:,:,ind),3)+(~hphie(:,:,i+1));
                    vlfe_phi(i,:,j) = sum(cond==(i-1))/ns;
                end
            end
            cond = sum(hphie(:,:,1:i),3);
            vlfe_phi(i,:,i+1) =-sum(cond==i)/ns;
        else
            cond = sum(hphie(:,:,1:i),3);
            vlfe(i,:) = sum(cond==i)/ns;
            for j=1:i
                ind = [setdiff(1:i,j)];
                cond = sum(hphie(:,:,ind),3);
                vlfe_phi(i,:,j) = sum(cond==(i-1))/ns;
            end
        end
        VC(i,opc) = sum(vlfe(i,:))/noe;
        if opc<OP.nRelax(i)
            Lambda(i)=OP.Mu(i)*(VC(i,opc)-VC(i,1)+(VC(i,1)-OP.VC0(i))*opc/OP.nRelax(i)); % Lambda calcultion in augmented Lagrngian
        else
            Lambda(i) =Lambda(i)+OP.Gamma(i)*(VC(i,opc)-OP.VC0(i));
            OP.Gamma(i) = min(OP.Gamma(i)+OP.dGamma(i),OP.maxGamma(i));
        end
    end
    ErModel = sum((vlfe.^MP.P).*MP.Emax')+(1-sum(vlfe))*MP.Emin;% ErModel=Ersatz material model and stiffness calculation
    % Assembling of stiffness matrix
    sK = reshape(Ke(:)*(ErModel.*ones(1,noe)),8*8*noe,1);
    K = sparse(iK,jK,sK);
    K = (K+K')/2;
    Dis(ukdis)=K(ukdis,ukdis)\F(ukdis);% Solving system of equlibrium equations
    ElemComp=sum(0.5*(Ke*Dis(EDofMat)').*(Dis(EDofMat)'.*ErModel));
    Compliance(opc) = sum(ElemComp);
    ErModel_phi = MP.P*sum(vlfe_phi.*(vlfe.^(MP.P-1)).*MP.Emax');
    VC_phi = sum(Lambda.*vlfe_phi);  
    for i=1:GP.noph
        ElemComp_phi=-sum(0.5*(Ke*Dis(EDofMat)').*(Dis(EDofMat)'.*ErModel_phi(:,:,i)));
        Comp_phi = sparse(ElNodes,ones(noe,4),0.25*ElemComp_phi'.*ones(noe,4));% Comp_phi = objective function sensitivity wrt Phi
        A_phi = sparse(ElNodes,ones(noe,4),0.25*VC_phi(:,:,i)'.*ones(noe,4));% A_phi = area function sensitivity wrt Phi
        V = Comp_phi/mean(abs(Comp_phi))+A_phi;% V = Boundary velocity in LS method
        T = (T1/OP.delta_T+OP.Tho*T2);
        Yy = (T1*(Phi(:,:,i)/OP.delta_T-V));
        Phi(:,:,i)=T\Yy;
        Phi(:,:,i) = min(max(Phi(:,:,i),-1),1);
        right_edge_nodes = (NodeCoord(:, 1) >= (GP.lx - GP.els)) & (NodeCoord(:, 1) <= GP.lx) &(NodeCoord(:, 2) >= (GP.ly-GP.cutout_size)/3) & (NodeCoord(:, 2) <= 2*(GP.ly-GP.cutout_size)/3);
        left_edge_nodes = (NodeCoord(:, 1) <=  GP.els) & (NodeCoord(:, 1) >= 0) &(NodeCoord(:, 2) >= 0) & (NodeCoord(:, 2) <= GP.ly);
        boundary_nodes1 = right_edge_nodes; 
        boundary_nodes2 = NodeCoord(:, 1) >= 0 & NodeCoord(:, 1) <= 24 & NodeCoord(:, 2) >= (GP.ly-GP.els) & (NodeCoord(:, 2) <= GP.ly);
        if i>1
            ind=Phi(:,:,i-1)<0;
            Phi(ind,:,i)=Phi(ind,:,i-1)-0.8;  % 减小偏移量
        end
       
        Phi(boundary_nodes1,:,1)=0.55;  % 减小初始值
        Phi(boundary_nodes1,:,2)=0.45;
        Phi(boundary_nodes1,:,3)=0.35;
        Phi(boundary_nodes1,:,4)=0.25;
        Phi(boundary_nodes1,:,5)=0.15;
        Phi(boundary_nodes1,:,6)=0.05;
        Phi(boundary_nodes1,:,7)=0.03;
        Phi(boundary_nodes1,:,8)=0.02;
        Phi(boundary_nodes1,:,9)=0.01;
        
        Phi(boundary_nodes2,:,1)=0.08;  % 调整第二个边界条件
        Phi(boundary_nodes2,:,2)=-0.08;
        Phi(boundary_nodes2,:,3)=-0.15;
        Phi(boundary_nodes2,:,4)=-0.22;
        Phi(boundary_nodes2,:,5)=-0.28;
        Phi(boundary_nodes2,:,6)=-0.35;
        Phi(boundary_nodes2,:,7)=-0.42;
        Phi(boundary_nodes2,:,8)=-0.45;
        Phi(boundary_nodes2,:,9)=-0.47;
        for j = 1:size(NodeCoord,1)% 遍历所有节点
            if NodeCoord(j,1) >= (GP.lx - GP.cutout_size) && NodeCoord(j,2) >= (GP.ly - GP.cutout_size)%判断节点x and y坐标是否在挖掉区域
                Phi(j,1,i) = -1;  % 将挖掉区域的Phi设为-1
            end
        end
    end
    
       % 定义材料颜色和名称
    materialColors = [
    0.000, 0.447, 0.698;   % blue (kept from original start)

    0.850, 0.325, 0.098;   % red-orange
    0.494, 0.184, 0.556;   % purple
    0.301, 0.745, 0.933;   % light blue
    0.466, 0.674, 0.188;   % green
    0.635, 0.078, 0.184;   % dark red
    0.600, 0.600, 0.600;   % gray
    0.098, 0.098, 0.439;   % dark slate blue

    0.000, 0.000, 0.000    % black (kept from original end)
];
    
    materialNames = { ...
    'Mo45Ta35V20', ...
    'Mo65Ta5V30', ...
    'W10Ta30Mo45V15', ...
    'W20Ta25Mo45V10', ...
    'W35Ta20Mo45', ...
    'W50Ta15Mo35', ...
    'W60Ta15Mo25', ...
    'W75Ta15Mo10', ...
    'W95Ta5' ...
};
    
    figure(97); clf
    if OP.PlotResult == 1
        hold on
        axis equal
        for i = 1:GP.noph
            x_fine = linspace(0, double(nelx+1), double(nelx+1));
            y_fine = linspace(0, double(nely+1), double(nely+1));
            [x, y] = meshgrid(x_fine, y_fine);
            z = reshape(Phi(:,:,i), (nelx+1), (nely+1))';
            contourf(x, y, z, [0, 0], 'FaceColor', materialColors(i,:));
        end
    
        % 添加 legend 替代 colorbar
        for i = 1:GP.noph
            h(i) = patch(NaN, NaN, materialColors(i,:));
        end
        legend(h, materialNames, 'Location', 'eastoutside');
    
        hold off
        drawnow
    end

    
    % --- Volume evolution plot ---
    figure(66); clf
    hold on
    iterations = 1:opc;
    for i = 1:GP.noph
        plot(iterations, VC(i, iterations), ...
             'LineWidth', 2, ...
             'Color', materialColors(i,:), ...
             'DisplayName', materialNames{i});
    end
    xlabel('Iteration', 'FontSize', 12);
    ylabel('Volume Fraction', 'FontSize', 12);
    title('Volume Evolution of Each Phase', 'FontSize', 14);
    legend('Location', 'best');
    grid on;
    box on;
    set(gca, 'FontSize', 12);
    hold off

    
    figure(68); clf;
    iterations = 1:opc;
    plot(iterations, Compliance(iterations), ...
     'LineWidth', 2, ...
     'Color', [0 0 0], ...  % Use a fixed color, e.g., black
     'DisplayName', 'Compliance');
    xlabel('Iteration', 'FontSize', 12);
    ylabel('Compliance', 'FontSize', 12);
    title('Compliance Evolution', 'FontSize', 14);
    legend('Location', 'best');
    grid on;
    box on;
    set(gca, 'FontSize', 12);
    total_cost = sum(VC(:, opc)' .* MP.cost);          % 每种材料的成本之和
    total_mass = sum(VC(:, opc)' .* MP.density);       % 每种材料的质量之和
    fprintf('opc: %d, Volumes: [%s], Compliance: %.4f\n', ...
    opc, strjoin(compose('%.4f', VC(:, opc)'), ', '), Compliance(opc));
    %fprintf('opc: %d\n', opc);
    %fprintf('Volumes[Al, Ag, Cu, TS]: %s\n', mat2str(VC(:, opc)', 4)); % 每种材料的体积
    %fprintf('Compliance: %.4f\n', Compliance(opc));
    %fprintf('Total Cost: %.2f\n', total_cost);
    %fprintf('Total Mass: %.2f\n', total_mass);
    if opc > 1
        rel_change = abs(Compliance(opc) - Compliance(opc-1))/Compliance(opc-1);
        if rel_change > 0.1
            OP.delta_T = max(OP.delta_T * 0.5, 0.001);  % 如果变化太大，减小时间步长
        elseif rel_change < 0.01
            OP.delta_T = min(OP.delta_T * 1.1, 0.05);  % 如果变化太小，适当增加时间步长
        end
        
        % 添加收敛检查
        %if opc > 10
            %recent_changes = abs(diff(Compliance(opc-9:opc)))./Compliance(opc-9:opc-1);
            %if mean(recent_changes) < 0.001
                %fprintf('Converged after %d iterations\n', opc);
                %break;
            %end
        %end
    end
    opc=opc+1;
end
total_time = toc;  % 结束计时
fprintf('Total execution time: %.2f seconds\n', total_time);
% close(v)
%% Mesh generation function
function [ndof,noe,NodeCoord,EDofMat,ElNodes,NodeRep,Phi,IC,nelx,nely] = MeshGeneration(lx,ly,els,noph,cutout_size)
    coordinates = cell(1,1);
    NodeCoord = [];
    nelx = int64(lx/els); 
    nely = int64(ly/els); 
    nonx = nelx+1; % number of nodes in x direction
    nony = nely+1; % number of nodes in y direction
    x_node = linspace(0,lx,nonx);
    y_node = linspace(0,ly,nony);
    x = repmat(x_node',nony,1);
    y = repmat(kron(y_node',ones(nonx,1)),1);
    coordinates = [x,y];
    NodeCoord = [NodeCoord;coordinates];%每个节点的x, y 坐标，维度为（nonx*nony,2),nonx*nony是节点数  
    NodeCoord = unique(NodeCoord,'stable','rows');
    
    % 创建L形设计域
    Phi = 0.1*ones(size(NodeCoord,1),1,noph);
    
    % 标记要挖掉的区域
    for i = 1:noph
        for j = 1:size(NodeCoord,1)% 遍历所有节点
            if NodeCoord(j,1) >= (lx - cutout_size) && NodeCoord(j,2) >= (ly - cutout_size)%判断节点x and y坐标是否在挖掉区域
                Phi(j,1,i) = -1;  % 将挖掉区域的Phi设为-1
            end
        end
    end
    
    ElNodes = [];     
    IC = cell(1,1);
    nodes = (1:size(coordinates,1))';
    [~,ic,iC] = intersect(coordinates,NodeCoord,'stable','rows');% ic 表示坐标在NodeCoord中的行索引，iC 表示NodeCoord中的行索引，都表示具体有多少个节点
    nodes(ic) = iC;IC=iC;
    a = double(repmat([0,1,nelx+1,nelx+2],nelx,1));
    b = double(repmat(a,nely,1)) + double(kron(double(0:nely-1)',ones(nelx,1,'double')));
    eleNode = double(b) + double(1:nelx*nely)';%每一个单元的节点索引，维度为（nelx*nely,4)
    eleNode = nodes(eleNode);
    ElNodes = [ElNodes;eleNode]; 
    EDofMat = kron(ElNodes,[2,2])+repmat([-1,0],1,4);% 每个单元的8个自由度，维度为（nelx*nely,8)
    noe = sum(nelx.*nely);
    nonodes = max(ElNodes,[],'all');%得到的节点索引的最大值，代表有多少个节点
    ndof = 2*nonodes;%每个节点有两个自由度，ndof即代表总共有多少个自由度
    NodeRep=groupcounts(ElNodes(:));% NodeRep = reppetition of each node in all element;
end
%% Boundary conditions implementation function
function [F,Dis,FixedNodes,ukdis,kdis] = BoundaryConditionsImplementation(ndof,NodeCoord)
    Dis = nan(ndof,1);  % 位移向量
    F = zeros(ndof,1);  % 力向量
    
    % 顶部边界节点（完全固定）
    top_nodes = NodeCoord(:,2) == max(NodeCoord(:,2));  % 顶部边界节点
    FixedNodes = find(top_nodes);     % 固定节点的编号
    Dis(FixedNodes*2-1) = 0;  % x方向位移为0
    Dis(FixedNodes*2) = 0;    % y方向位移为0
    
    % 右边界中点施加向下的力（在0和(ly - cutout_size)之间）
    right_nodes = NodeCoord(:,1) == max(NodeCoord(:,1));  % 右边界节点
    mid_y = (60 - 36)/2;  % 在0和(ly - cutout_size)之间的中点
    load_nodes = find(right_nodes & abs(NodeCoord(:,2) - mid_y) < 1);  % 中点附近的节点
    F(load_nodes*2) = -100;    % y方向力为-1
    ukdis = isnan(Dis);  % 未知位移的索引
    kdis = ~isnan(Dis);  % 已知位移的索引
end
%% Plane element stiffness calculation function
function [Ke,Vole,Bmat,Emat] = PlaneElementStiffnessCalculation(ngp,els,v)
    eXcor = [0,els,0,els]';
    eYcor = [0,0,els,els]';
    Emat = 1/(1-v^2)*[1,v,0;v,1,0;0,0,(1-v)/2];% Elasticity matrix by cosidering 1 as modulus of elasticity
    Bmat = zeros(3,8);
    Volg = zeros(ngp);% Volume at gussian points
    Vole = 0;% Vole = Volume for each element
    Ke = zeros(8);% Ke = stiffness matrix for element
    [gp,wgp]=makegaussianpoint(ngp);% gauss points and their wheights
    for j=1:ngp
        s = gp(j);
        for k=1:ngp
            r = gp(k);
            N_r = 0.25*[-(1-s),(1-s),-(1+s),(1+s)];
            N_s = 0.25*[-(1-r),-(1+r),(1-r),(1+r)];
            X_r = N_r*eXcor; Y_r = N_r*eYcor;% Calcualtion derivitave of global coordiantes wrt local coordiantes
            X_s = N_s*eXcor; Y_s = N_s*eYcor;
            J = [X_r,Y_r;X_s,Y_s];% Jacobian matrix calculation
            N_X_Y = J\[N_r;N_s];
            Bmat(1,1:2:end)=N_X_Y(1,:);
            Bmat(2,2:2:end)=N_X_Y(2,:);% For Normal strain at y direction
            Bmat(3,1:2:end)=N_X_Y(2,:);% For Shear strain at x and y direction
            Bmat(3,2:2:end)=N_X_Y(1,:);
            Volg(k,j) = det(J)*wgp(k)*wgp(j);% Volume at gussian points
            Vole = Vole+Volg(k,j);% Volume for each element = sum(Arg)
            Ke = Ke+Bmat'*Emat*Bmat*Volg(k,j);% Stiffness for each element = sum(Stiffness at each gussian point)
        end
    end
end
%% T1 and T2 Generator function
function [T1e,T2e] = RD_T1T2Generator(ngp,els)
    eXcor = [0,els,0,els]';
    eYcor = [0,0,els,els]';
    Volg = zeros(ngp,ngp);% Volume at gussian points
    T1e = zeros(4);
    T2e = zeros(4);
    [gp,wgp]=makegaussianpoint(ngp);
    for j=1:ngp
        s = gp(j);
        for k=1:ngp
            r = gp(k);
            N=0.25*[(1-r)*(1-s),(1+r)*(1-s),(1-r)*(1+s),(1+r)*(1+s)];
            N_r = 0.25*[-(1-s),(1-s),-(1+s),(1+s)];
            N_s = 0.25*[-(1-r),-(1+r),(1-r),(1+r)];
            X_r = N_r*eXcor; Y_r = N_r*eYcor;
            X_s = N_s*eXcor; Y_s = N_s*eYcor;
            J = [X_r,Y_r;X_s,Y_s];
            N_X_Y = J\[N_r;N_s];
            Volg(k,j) = det(J)*wgp(k)*wgp(j);
            T1e = T1e+N'*N*Volg(k,j);% T1 and T2 for each element = sum(Stiffness at each gussian point)
            T2e = T2e+N_X_Y'*N_X_Y*Volg(k,j);
        end
    end
end
function [r,W]=makegaussianpoint(ng)
    switch ng
            case 1
                    r=0;W=2;
            case 2
                    r=[-0.577,0.577];W=[1,1];
            case 3
                    r=[-0.774 0 0.774];W=[0.555 0.888 0.555];
            case 4
                    r=[-0.861 -0.34 0.34 0.861];W=[0.348 0.652 0.652 0.348];
            case 5
                    r=[-0.906 -0.538 0 0.538 0.906];W=[0.237 0.478 0.568 0.478 0.237];
    end
end