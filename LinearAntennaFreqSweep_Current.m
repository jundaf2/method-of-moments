%% plot the current distribution at 1:5 *Frequency
clear all;close all;
Frequency=1e9;
freq_point=5;
Ny = 100; % number of segments + 1
Vm = 1; % delta-gap voltage
%% Basic parameters
c = 299792458;                       
Mu = 4.0e-7*pi;                      
Epsilon = 1.0/(Mu*c*c);              
NumGauss=5; % number of Gauss Quadrature used

%% Input parameters
lambda = 1/Frequency*c;
w = 2.0*pi*Frequency;                % angular freq
k = w*sqrt(Mu*Epsilon);              % wave number
Eta = sqrt(Mu/Epsilon);              % wave impedance, aka Z0
% coordinates of the corner points 
p0=[0 0];
p1=[0.5*lambda 0];
WireRadius=1e-3*lambda;  % radius of thin wire
line_length=sqrt((p1(2)-p0(2))^2+(p1(1)-p0(1))^2);

%% Freq sweep
Frequency_List= (1:freq_point)*Frequency; % 1f0 2f0 3f0 4f0
Zin_NumericalResultList=zeros(1,freq_point);
w_List = 2.0*pi*Frequency_List;                % angular freq
k_List = w_List*sqrt(Mu*Epsilon);              % wave number
%% Generate the mesh
ty = linspace(0,1,Ny); % equally spaced segments
% line antenna geometry
x1 = p0(1)+(p1(1)-p0(1))*ty;
y1 = p0(2)+(p1(2)-p0(2))*ty;

TotalSegments=(Ny-1); % N in Prof. Jin's book
Delta_Gap_Length = line_length/TotalSegments;
TotalElements = TotalSegments-1; 

for i=1:Ny
    NodeX(i)=x1(i);
    NodeY(i)=y1(i);
end

% index of segment
first=1;
mid = floor(TotalSegments/2)+1;
final=TotalElements;

SourceId = mid;

%% Node List counted by element, also contains segment information --> 1:2 (+), 2:3(-)
ElementNodeList = zeros(TotalElements, 3);
for i=1:TotalElements
    ElementNodeList(i,1) = i;
    ElementNodeList(i,2) = i+1;
    ElementNodeList(i,3) = i+2;
end
    
figure(1);hold on;
plot(NodeX,NodeY,'k*-','linewidth',2);title('thin wire geometry');

plot(NodeX(SourceId:SourceId+1),NodeY(SourceId:SourceId+1),'bo-','linewidth',2);

legend(['Segment'],['Source']);
hold off;
%% basis functions and gauss integral
[GaussL, GaussW] = lgwt(NumGauss, 0.0, 1);  % Gauss integration weights and locations (Gauss-Legendre)

tri(1,:) = 1 - GaussL;      % ascending triangle function for numerical integration 
tri(2,:) = GaussL;          % descending triangle function for numerical integration

nds = 20;                   % number of numerical integration points for self terms
%% construct the Z^{EJ}_{mn} Matrix
A = zeros(TotalElements,TotalElements,freq_point);   % zero matrix Ab=rhs

tic
% loop all test functions
for m_tri = 1:TotalElements

    ms(1) = ElementNodeList(m_tri,1);  % m tri start node
    ms(2) = ElementNodeList(m_tri,2);  % m tri peak node
    ms(3) = ElementNodeList(m_tri,3);  % m tri end node

    % node coordinate of m tri
    xms(1) = NodeX(ms(1)); 
    yms(1) = NodeY(ms(1));
    xms(2) = NodeX(ms(2)); 
    yms(2) = NodeY(ms(2));
    xms(3) = NodeX(ms(3)); 
    yms(3) = NodeY(ms(3));
    DeltaLM(1)=sqrt((xms(1)-xms(2))^2+(yms(1)-yms(2))^2);
    DeltaLM(2)=sqrt((xms(3)-xms(2))^2+(yms(3)-yms(2))^2);

    ds(1) = DeltaLM(1)/nds;
    ds(2) = DeltaLM(2)/nds;
    sx(1,:) = linspace(0.5*ds(1), DeltaLM(1) - 0.5*ds(1), nds); % approximate triangular basis function as a constant
    sx(2,:) = linspace(0.5*ds(2), DeltaLM(2) - 0.5*ds(2), nds);           
    tris(1,:) = sx(1,:) /DeltaLM(1);              % ascending triangle function for numerical integration in self term
    tris(2,:) = (DeltaLM(2) - sx(2,:))/DeltaLM(2);   % descending triangle function for numerical integration in self term

    % loop all "source" triangles
    for n_tri = 1:TotalElements % number of triangular basis function
        
        ns(1) = ElementNodeList(n_tri,1);  % n tri start node
        ns(2) = ElementNodeList(n_tri,2);  % n tri peak node
        ns(3) = ElementNodeList(n_tri,3);  % n tri end node
        % node coordinate of n tri
        xns(1) = NodeX(ns(1));
        yns(1) = NodeY(ns(1));
        xns(2) = NodeX(ns(2));
        yns(2) = NodeY(ns(2));
        xns(3) = NodeX(ns(3));
        yns(3) = NodeY(ns(3));
        DeltaLN(1)=sqrt((xns(1)-xns(2))^2+(yns(1)-yns(2))^2);
        DeltaLN(2)=sqrt((xns(3)-xns(2))^2+(yns(3)-yns(2))^2);
        
        for im = 1:2    % loop testing segments
            % positive or negative
            if im == 1 %first segment +
                sm = 1.0;
            else
                sm = -1.0;%second segment -
            end
            for in=1:2 % loop source segments
                % positive or negative derivative
                if in == 1
                    sn = 1.0;
                else
                    sn = -1.0;
                end
                if ms(im) == ns(in)% self-term
                    % equ (10.3.41)
                    % mid point approximation for the triangular function
                    DeltaL = DeltaLN(im);
                    Phi = (DeltaL/2/pi*(log(2*DeltaL/WireRadius)-1)-1j*k_List/4/pi*DeltaL^2);
                    s1 = tris(in,nds/2)*tris(im,nds/2).*Phi*1j.*k_List*Eta;
                    s2 = -sn*sm.*Phi/(DeltaL^2).*1j*Eta./k_List;
                    A(m_tri,n_tri,:) = A(m_tri,n_tri,:) + reshape(s1+s2,[1 1 freq_point]);
                    
                else % non-self term
                    sv = 0.0;
                    ss = 0.0;
                    SegmentNVecM = [ xms(im+1)-xms(im), yms(im+1)-yms(im) ];
                    SegmentNVecM = SegmentNVecM/norm(SegmentNVecM); % normalized

                    SegmentNVecN = [ xns(im+1)-xns(im), yns(im+1)-yns(im) ];
                    SegmentNVecN = SegmentNVecN/norm(SegmentNVecN); % normalized

                    xm = xms(im) + DeltaLM(im)*(1.0-GaussL)*SegmentNVecM(1);  % gauss quadrature points on m triangle 
                    ym = yms(im) + DeltaLM(im)*(1.0-GaussL)*SegmentNVecM(2);

                    xn = xns(in) + DeltaLN(in)*(1.0-GaussL)*SegmentNVecN(1); % gauss quadrature points on n triangle 
                    yn = yns(in) + DeltaLN(in)*(1.0-GaussL)*SegmentNVecN(2);

                    for i_gm = 1:NumGauss 
                        for i_gn = 1:NumGauss
                            r = sqrt((xm(i_gm) - xn(i_gn))^2 + (ym(i_gm) - yn(i_gn))^2 + WireRadius^2);
                            expterm = exp(-1j*reshape(k_List,[1 1 freq_point]).*repmat(r,1,1,freq_point))./r/4/pi; % 3d-freespace green function
                            sv = sv + DeltaLM(im)*DeltaLN(in)*GaussW(i_gm)*GaussW(i_gn)*tri(im,i_gm)*tri(in,i_gn)*expterm*(SegmentNVecM*SegmentNVecN'); % vector f_m and vector f_n 
                            ss = ss + GaussW(i_gm)*GaussW(i_gn)*expterm;
                        end
                    end
                    A(m_tri,n_tri,:) = A(m_tri,n_tri,:) + reshape(k_List.*k_List.*reshape(sv,[1 freq_point 1]) - sm*sn.*reshape(ss,[1 freq_point 1]),[1 1 freq_point]).*reshape(1j*Eta./k_List,[1 1 freq_point]);
                end
            end
        end
    end

end
%% excitation in delta-gap form (compute right-hand side )
% rhs = zeros(TotalElements,1); 
rhs = zeros(TotalElements,freq_point); 
for num_tri = 1:TotalElements
    s(1) = num_tri;
    s(2) = num_tri + 1;
    for im=1:2 
        if s(im) == SourceId
            rhs(num_tri,:)=Vm/2;
        end
    end
end
toc;
%% build current distribution from solved parameters
CurrentList = cell(1,freq_point);
CurrentInterpList = cell(1,freq_point);
for i=1:freq_point
    [L,U]=lu(A(:,:,i));
    b=U\(L\rhs(:,i));%A(:,:,i)\rhs(:,i);%
    % compute input impedance at the middle segment
    Zin = Vm / (0.5*(b(SourceId) + b(SourceId+1)));
    Zin_NumericalResultList(i)=Zin;
    
    Current = zeros(TotalSegments,1); % using the average of neiboring two basis function
    for j=1:TotalSegments
        if j == 1
            Current(j) = 0.5*b(j);
        elseif j == TotalSegments
            Current(j) = 0.5*b(TotalSegments - 1);
        else
            Current(j) = 0.5*(b(j-1) + b(j));
        end
    end
    CurrentList{i}=Current;
    
    CurrentInterp = zeros(TotalSegments*NumGauss,1); % interpolation using basis function
    for num_tri = 1:TotalElements
        s(1) = num_tri-1; % segment id - 1 
        s(2) = num_tri; % segment id - 1 
        for im=1:2 
            for i_gm = 1:NumGauss 
                CurrentInterp(s(im)*NumGauss+i_gm)=CurrentInterp(s(im)*NumGauss+i_gm)+tri(im,i_gm)*b(num_tri);
            end
        end
    end
    CurrentInterpList{i}=CurrentInterp;
end

%%
figure(2);
subplot(121)
hold on;
plot(line_length*(1:TotalSegments*NumGauss)/(TotalSegments*NumGauss),abs(CurrentInterpList{1}),'-r','Linewidth',3);
plot(line_length*(1:TotalSegments*NumGauss)/(TotalSegments*NumGauss),abs(CurrentInterpList{2}),'-g','Linewidth',3);
plot(line_length*(1:TotalSegments*NumGauss)/(TotalSegments*NumGauss),abs(CurrentInterpList{3}),'-b','Linewidth',3);
plot(line_length*(1:TotalSegments*NumGauss)/(TotalSegments*NumGauss),abs(CurrentInterpList{4}),'-c','Linewidth',3);
plot(line_length*(1:TotalSegments*NumGauss)/(TotalSegments*NumGauss),abs(CurrentInterpList{5}),'-m','Linewidth',3);
grid on;xlabel('physical length');ylabel('Magnitude');
legend('f_0','2f_0','3f_0','4f_0','5f_0');
set(gca,'fontsize',24);
hold off
subplot(122)
hold on;
plot(line_length*(1:TotalSegments*NumGauss)/(TotalSegments*NumGauss),angle(CurrentInterpList{1})/pi*180,'-r','Linewidth',3);
plot(line_length*(1:TotalSegments*NumGauss)/(TotalSegments*NumGauss),angle(CurrentInterpList{2})/pi*180,'-g','Linewidth',3);
plot(line_length*(1:TotalSegments*NumGauss)/(TotalSegments*NumGauss),angle(CurrentInterpList{3})/pi*180,'-b','Linewidth',3);
plot(line_length*(1:TotalSegments*NumGauss)/(TotalSegments*NumGauss),angle(CurrentInterpList{4})/pi*180,'-c','Linewidth',3);
plot(line_length*(1:TotalSegments*NumGauss)/(TotalSegments*NumGauss),angle(CurrentInterpList{5})/pi*180,'-m','Linewidth',3);
grid on;xlabel('physical length');ylabel('Phase Angle');
legend('f_0','2f_0','3f_0','4f_0','5f_0');
set(gca,'fontsize',24);
hold off

figure(4);
subplot(121)
hold on;
plot(line_length*(1:TotalSegments*NumGauss)/(TotalSegments*NumGauss),abs(CurrentInterpList{1}),'-r','Linewidth',4);
grid on;xlabel('physical length');ylabel('Magnitude');
set(gca,'fontsize',24);
hold off
subplot(122)
hold on;
plot(line_length*(1:TotalSegments*NumGauss)/(TotalSegments*NumGauss),angle(CurrentInterpList{1})/pi*180,'-b','Linewidth',4);
grid on;xlabel('physical length');ylabel('Phase Angle');
set(gca,'fontsize',24);
hold off

figure(3);hold on;
AbsReCur = abs(real(CurrentList{1}));
IntAbsReCur = AbsReCur ./ (double(vpa(max(AbsReCur),2))/(1+floor(max(AbsReCur) ./ 10.^floor(log10(max(AbsReCur))))));
max_color_value = 100*ceil(max(IntAbsReCur));
jet_color = colormap(jet(ceil(max_color_value)));
title('Real Part of Current Along Wire');
for i=1:TotalSegments-1
        selected_color = jet_color(ceil(100*IntAbsReCur(i)),:);
        plot(NodeX(i:i+1),NodeY(i:i+1),'*-','linewidth',2,'color',selected_color);
end
colorbar;
hold off;