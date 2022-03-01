%% sweep for Zin * 50 size 
clear all;close all;
Frequency=1e9;
segment_point_fac10 = 60;
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
Eta = sqrt(Mu/Epsilon);              % wave impedance
% coordinates of the corner points 
p0=[0 0];
p1=[0.5*lambda 0];
WireRadius=1e-3*lambda;  % radius of thin wire
line_length=sqrt((p1(2)-p0(2))^2+(p1(1)-p0(1))^2);
Ny_List=(1:segment_point_fac10)*10;
Zin_NumericalResultList=zeros(1,segment_point_fac10);
for seg_id=1:length(Ny_List)
    Ny=Ny_List(seg_id); % loop through different number of segments
    %% Generate the mesh
    ty = linspace(0,1,Ny) ;
    % line antenna geometry
    x1 = p0(1)+(p1(1)-p0(1))*ty;
    y1 = p0(2)+(p1(2)-p0(2))*ty;

    TotalSegments=(Ny-1);
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
    ElementNodeList(1,1) = 1;
    ElementNodeList(1,2) = 2;
    ElementNodeList(1,3) = 3;
    for i=2:TotalElements-2
        ElementNodeList(i,1) = i;
        ElementNodeList(i,2) = i+1;
        ElementNodeList(i,3) = i+2;
    end


    ElementNodeList(TotalElements-1,1) = TotalElements-1;
    ElementNodeList(TotalElements-1,2) = TotalElements;
    ElementNodeList(TotalElements-1,3) = TotalElements+1;

    ElementNodeList(TotalElements,1) = TotalElements;
    ElementNodeList(TotalElements,2) = TotalElements+1;
    ElementNodeList(TotalElements,3) = TotalElements+2;

    % middle of a segment, no use, unless need visualization
    MidNodeX = zeros(1,TotalSegments);
    MidNodeY = zeros(1,TotalSegments);
    for i = 1:TotalSegments-1
        MidNodeX(i) = 0.5*(NodeX(i)+NodeX(i+1)); 
        MidNodeY(i) = 0.5*(NodeY(i)+NodeY(i+1)); 
    end

    MidNodeX(TotalSegments) = 0.5*(NodeX(TotalSegments)+NodeX(1)); 
    MidNodeY(TotalSegments) = 0.5*(NodeY(TotalSegments)+NodeY(1)); 

    %% basis functions and gauss integral
    [GaussL, GaussW] = lgwt(NumGauss, 0.0, 1);  % Gauss integration weights and locations (Gauss-Legendre)

    tri(1,:) = 1 - GaussL;      % ascending triangle function for numerical integration 
    tri(2,:) = GaussL;          % descending triangle function for numerical integration

    nds = 20;                   % number of numerical integration points for self terms
    %% construct the Z^{EJ}_{mn} Matrix
    A = zeros(TotalElements);   % zero matrix Ab=rhs

    % loop all test functions
    for m_tri = 1:TotalElements % 试函数个数 Gibson的N

    ms(1) = ElementNodeList(m_tri,1);  % m tri start
    ms(2) = ElementNodeList(m_tri,2);  % m tri peak
    ms(3) = ElementNodeList(m_tri,3);  % m tri end

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
    sx(1,:) = linspace(0.5*ds(1), DeltaLM(1) - 0.5*ds(1), nds);    
    sx(2,:) = linspace(0.5*ds(2), DeltaLM(2) - 0.5*ds(2), nds);           
    tris(1,:) = sx(1,:) /DeltaLM(1);              % ascending triangle function for numerical integration in self term
    tris(2,:) = (DeltaLM(2) - sx(2,:))/DeltaLM(2);   % descending triangle function for numerical integration in self term

    % loop all "source" triangles
        for n_tri = 1:TotalElements
            ns(1) = ElementNodeList(n_tri,1);  % n tri start
            ns(2) = ElementNodeList(n_tri,2);  % n tri peak
            ns(3) = ElementNodeList(n_tri,3);  % n tri end
            % node coordinate of n tri
            xns(1) = NodeX(ns(1));
            yns(1) = NodeY(ns(1));
            xns(2) = NodeX(ns(2));
            yns(2) = NodeY(ns(2));
            xns(3) = NodeX(ns(3));
            yns(3) = NodeY(ns(3));
            DeltaLN(1)=sqrt((xns(1)-xns(2))^2+(yns(1)-yns(2))^2);
            DeltaLN(2)=sqrt((xns(3)-xns(2))^2+(yns(3)-yns(2))^2);
            DeltaL = DeltaLN(1);
            for im = 1:2    % loop testing segments

                % positive or negative derivative for scalar term
                if im == 1 % +
                    sm = 1.0;
                else
                    sm = -1.0;% -
                end

                for in=1:2 % loop source segments

                    % positive or negative derivative for 5.74 and 5.81 scalar term
                    if in == 1
                        sn = 1.0;
                    else
                        sn = -1.0;
                    end
                    if ms(im) == ns(in)% self-term
                        % equ (10.3.41)
                        % mid point approximation for the triangular function
                        DeltaL = DeltaLN(im);
                        Phi = (DeltaL/2/pi*(log(2*DeltaL/WireRadius)-1)-1j*k/4/pi*DeltaL^2);
                        s1 = tris(in,nds/2)*tris(im,nds/2).*Phi*1j.*k*Eta;
                        s2 = -sn*sm.*Phi/(DeltaL^2).*1j*Eta./k;
                        A(m_tri,n_tri) = A(m_tri,n_tri) + s1+s2;
                    else % non-self term
                        sv = 0.0;
                        ss = 0.0;
                        % equ 5.75

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
                                expterm = exp(-1i*k*r)/r/4/pi; % 3d-freespace green function
                                sv = sv + DeltaLM(im)*DeltaLN(in)*GaussW(i_gm)*GaussW(i_gn)*tri(im,i_gm)*tri(in,i_gn)*expterm*(SegmentNVecM*SegmentNVecN'); % vector f_m and vector f_n 
                                ss = ss + GaussW(i_gm)*GaussW(i_gn)*expterm;
                            end
                        end
                        A(m_tri,n_tri) = A(m_tri,n_tri) + (k*k*sv - sm*sn*ss)*1j*Eta/k;
                    end
                end
            end
        end

    end
    
    %% excitation in delta-gap form (compute right-hand side )
    rhs = zeros(TotalElements,1); 
    for num_tri = 1:TotalElements
        s(1) = num_tri;
        s(2) = num_tri + 1;
        for im=1:2 
            if s(im) == SourceId
                rhs(num_tri,:)=Vm/2;
            end
        end
    end
    %% build current distribution from solved parameters
    b=A\rhs; % [A]^(-1) * [rhs] 

    % compute input impedance at 1st segment
    Zin = Vm / (0.5*(b(SourceId) + b(SourceId+1)));
    Zin_NumericalResultList(seg_id)=Zin;
end
%%

figure(4);hold on;
plot(Ny_List,real(Zin_NumericalResultList),'-r','Linewidth',4);
plot(Ny_List,imag(Zin_NumericalResultList),'-b','Linewidth',4);
grid on;xlabel('Number of Segments');
legend('Resistance=Real[Z_{in}]','Reactance =Imag[Z_{in}]');
set(gca,'fontsize',24);
hold off