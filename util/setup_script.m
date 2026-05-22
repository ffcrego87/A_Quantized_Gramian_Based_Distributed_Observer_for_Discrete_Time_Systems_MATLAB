%% System
disp('Setup script')

% system parameters
nx     = NAg*dim;    % Total dimension of the system
ny     = dim;        % Dimension of local outputs

% Generate network
gen_networks

% Process and measurement noise covariances
W = eye(nx);
V = eye((2*NAg-1)*ny);

%% Global Control Cost

R = eye(NAg);
Q = eye(nx);

%% system matrices

%% A
A = cell(ifinal,1);
for i=1:ifinal
    Aloc = eye(dim)+Amag*randn(dim,dim);
    A{i} = kron(eye(NAg),Aloc); % Global dynamic matrix
end

%% C
Cglob = cell(ifinal,1);
C = cell(NAg,1);
for t=1:ifinal
    Cloc = eye(dim)+Cmag*randn(dim,dim);
    G=zeros(2,NAg);
    G(1,1) = 1;
    G(2,1) = 1;
    G(1,2) = -1;
    G(2,NAg) = -1;
    C{1} = kron(G,Cloc);
    for i=2:(NAg-1)
        G=zeros(2,NAg);
        G(1,i) = 1;
        G(2,i) = 1;
        G(1,i+1) = -1;
        G(2,i-1) = -1;
        C{i} = kron(G,Cloc);
    end
    G=zeros(1,NAg);
    G(NAg) = 1;
    C{NAg} = kron(G,Cloc);
    Cglob{t} = cell2mat(C); % Global output matrix
end

% yidxs
yidxs = cell(NAg,1);
yidxs{1} = 1:(2*ny);
for i=2:(NAg-1)
    yidxs{i} = (i-1)*2*ny+(1:(2*ny));
end
yidxs{NAg} = (NAg-1)*2*ny+(1:ny);

%% Tau
Tau=zeros(NAg,NAg);
Piprev = eye(NAg);
t=1;
while ~min(min(Piprev))
    Pinext = sign(Pi*Piprev);
    Tau=Tau+t*(Pinext-Piprev);
    t=t+1;
    Piprev=Pinext;
end

%% N
n_mat = zeros(NAg,NAg);
for i=1:NAg
    for j=1:NAg
        if i==j
            n_mat(i,j) = i;
        else
            idxs = find(Tau(:,j)==(Tau(i,j)-1));
            n_mat(i,j) = idxs(find(Pi(i,idxs)>0,1));
        end
    end
end

%% D
D_cell = cell(NAg,NAg);
for i=1:NAg
    for j=1:NAg
        auxdiag=zeros(size(Cglob{1},1),1);
        neighs_active = find(n_mat(i,:)==j);
        for k = neighs_active
            auxdiag(yidxs{k})=ones(length(yidxs{k}),1);
        end
        D_cell{i,j}=diag(auxdiag);
    end
end

%% E
E_cell = cell(NAg,1);
for i=1:NAg
    E_cell{i}=zeros(size(Cglob{1},1)*(kfix+1),size(Cglob{1},1));
    for k=0:kfix
        auxdiag = zeros(size(Cglob{1},1),1);
        ag_toi = find(k==Tau(i,:));
        for ag = ag_toi
            auxdiag(yidxs{ag})=ones(length(yidxs{ag}),1);
        end
        E_cell{i}(size(Cglob{1},1)*k+(1:size(Cglob{1},1)),:)=diag(auxdiag);
    end
end

%% Delta
Delta_cell = cell(NAg,1);
for i=1:NAg
    auxdiag = zeros(size(Cglob{1},1)*(kfix+1),1);
    for k=0:kfix
        Pik=Pi^k;
        for j=1:NAg
            auxdiag(size(Cglob{1},1)*k+yidxs{j},1)=beta^(k+1)*Pik(i,j)*ones(length(yidxs{j}),1);
        end
    end
    Delta_cell{i}=diag(auxdiag);
end

if check_equivalence
    %% Omegas
    Omegabar = cell(NAg,ifinal);
    Omegabarprev = cell(NAg,1);
    Chi = alpha*eye(nx);
    for t=2:ifinal
        Chi = beta*(A{t}^(-1))'*Chi*A{t}^(-1);
        for i = 1:NAg
            Omegabarprev{i}=zeros(nx,nx);
            Omegabar{i,t}=zeros(nx,nx);
        end
        for tau=kfix-(0:kfix)
            tc = t-tau-1;
            if tc>=2
                for i = 1:NAg
                    Ccur = Cglob{tc}(yidxs{i},:);
                    Omegabar{i,t}=Ccur'*Ccur;
                    for j = 1:NAg
                        Omegabar{i,t}=Omegabar{i,t}+Pi(i,j)*Omegabarprev{j};
                    end
                    Omegabar{i,t}=beta*(A{tc}^(-1))'*Omegabar{i,t}*A{tc}^(-1);
                end
                Omegabarprev=Omegabar(:,t);
            end
        end
        
        for i=1:NAg
            Omegabar{i,t}=Omegabar{i,t}+Chi;
        end
    end
    
    %% Omegabar Gramian
    OmegabarReal = cell(NAg,ifinal);
    T_cell = cell(NAg,1);
    for i=1:NAg
        T_cell{i} = zeros(size(Cglob{1},1),nx);
    end
    O_cell = cell(NAg,1);
    for i=1:NAg
        O_cell{i} = zeros(size(Cglob{1},1)*(kfix+1),nx);
    end
    Chi = alpha*eye(nx);
    
    for t=2:ifinal
        Chi = beta*(A{t}^(-1))'*Chi*A{t}^(-1);
        T_prev_cell = T_cell;
        O_prev_cell = O_cell;
        for i = 1:NAg
            % T update
            T_cell{i}=0*T_cell{i};
            for j=1:NAg
                if i==j
                    T_cell{i}=T_cell{i}+D_cell{i,j}*Cglob{t}*A{t}^(-1);
                else
                    T_cell{i}=T_cell{i}+D_cell{i,j}*T_prev_cell{j}*A{t}^(-1);
                end
            end
            
            % O update
            O_update = [zeros(size(Cglob{t}));O_prev_cell{i}]*A{t}^(-1);
            O_cell{i} = O_update(1:(size(Cglob{1},1)*(kfix+1)),:)+E_cell{i}*T_cell{i};
            
            OmegabarReal{i,t}=O_prev_cell{i}'*Delta_cell{i}*O_prev_cell{i}+Chi;
            if max(max(abs(OmegabarReal{i,t}-Omegabar{i,t})))>10^(-5)
                fprintf('disagreement at iteration %i, agent %i\n',t,i)
            end
        end
    end
end

%% Quantization setup
% These parameters are used by sim_sample.m to build adaptive quantization
% ranges and step sizes at each iteration.
quantization_setup = struct;
quantization_setup.use_quantization = use_quantization;
quantization_setup.nbits = nbits;
quantization_setup.gammaOmega = gammaOmega;
quantization_setup.gammaMsg = gammaMsg;
quantization_setup.eps_quant = eps_quant;