clc
clear all
 
%% Data
E=200e9;            %Youngs Modulus
mew=0.3;            %Poisson Ratio
h=.02;              %thickness
ndof=2;             %degree of freedom per node
endcoo=[40 0;60 0;60 10;40 10 ];  %End coordinates of each node
nnx=10; % No of nodes in X direction
nny=5; % No of nodes in Y direction
ngp=3;  % No of Gauss Quadrature Points




%% Initialization
nnpe=4; %no of nodes per element
lx=endcoo(3,1)-endcoo(1,1);
ly=endcoo(3,2)-endcoo(1,2);
lex=lx/(nnx-1);         %Length of element in x dir
ley=ly/(nny-1);         %Length of element in x dir
tnn=nnx*nny;            %total no of nodes    
nex=nnx-1;              %No of elements in X direction
ney=nny-1;              %No of elements in y direction
tdof=tnn*ndof;          %total degree of freedom = number of nodes*dof per node
ne=2*2*(nex+ney); % Total No of elements
fg=zeros(tdof,1);       %force matrix initialization
kg=zeros(tdof,tdof);    %Global Elemental stiffness matrix generation
%calling function gaussquads for the weights and abssissa value. 
[weights,quadpts]=gaussquads(ngp);
weights=double(weights); %since the function used syms
quadpts=double(quadpts); %Zeta and Eta values
%D=(E/(1-mew^2))*[1 mew 0;mew 1 0; 0 0 (1-mew)/2];
D=(E/((1+mew)*(1-(2*mew))))*[1-mew mew 0 mew;mew 1-mew 0 mew; 0 0 (1-2*mew)/2 0;mew mew 0 1-mew];

%% Force Input and Boundary Conditions
%Force inout
%fg(49)=10000;fg(39)=20000;fg(29)=20000;fg(19)=20000;fg(9)=10000; %Force input in nodes
%Input force corresponding to coordinate system. 
%eg., Downward force is taken negative


%bodyforce=[100;200]; %bodyforce [bx;by];
 tractionx=2;tractiony=0; %traction force

 traction=[1,3];
 tn=zeros(1,ne);
 k=1;
 for i=traction
     tn(i)=traction(k);
     k=k+1;
 end
%% Boundary Condition
%bc(1)=;... % Manual data input
bc1=1:(nnx):tnn;%Boundary Conditions fixed nodes
for i=1:nny*size(bc1)
    bc(2*i-1)=2*bc1(i)-1;
    bc(2*i)=2*bc1(i);
end

%% Nodal Positioning
Q = 1:tnn;                      %Numbering elements
Q=reshape( Q,[ nnx nny]);       %Sorting element numbers 
Q=flip(Q');
Qx1=endcoo(1,1):lex:endcoo(2,1); 
Qy1=endcoo(1,2):ley:endcoo(3,2);
Qx=repmat(Qx1,[1, nny]);        %Qx- X position of nodes     
Qy=[];                          %Qy- Y position of nodes
for i=1:nny
a=repmat(Qy1(i) , [1,nnx]);
Qy=[Qy , a];
end
%  Qx=[40 50 60 40 50 60 40 50 60]; %X Coordinates of nodes 
%  Qy=[0 0 0 5 5 5 10 10 10]; %Y Coordinates of nodes 

%%  Quadrilateral Mesh generation
%Creating Quadrilateral connectivity
nnx1=nnx;
if nny>nnx
    nnx1=nny;
end
k=nnx:nnx:(nny*(nnx1-2));
quads=zeros(ne,nnpe);
for i=1:(tnn-(nnx+1))
    for j=1
        if(i~=k)
        quads(i,j)=i;
        quads(i,j+1)=i+1;
        quads(i,j+2)=i+nnx+1;
        quads(i,j+3)=i+nnx;
        end
    end
end
quads( all(~quads,2), : ) = [];

 %% Connection Matrix Generation
 for i=1:size(quads) 
       for j=1:nnpe %No of nodes per element
        n=quads(i,j);
        m=n*2-1;
        conn(i,2*j-1)=m;
        conn(i,2*j)=m+1; 
       end
 end
 
%% Computation
for i=1:size(quads)
     ke(i).element=zeros((nnpe*ndof),(nnpe*ndof));
     for j=1:nnpe
         X(j)=Qx(quads(i,j));
         Y(j)=Qy(quads(i,j));
     end
     xlim([0 2*max(Qx)]);
     ylim([0 2*max(Qy)]);
    
    plot(X,Y)%plots the graph for checking the connectivity
    hold on
     r1=X(1);r2=X(4);% since acting left side
     a=(2*r1+r2)/6;b=(r1+2*r2)/6; 
     
 if i==tn(i);
        fet(i).element=2*pi*ley*[a*tractionx a*tractiony 0 0 0 0 b*tractionx b*tractiony]';
 
 else
      fet(i).element=zeros(8,1);
 end

    for j=1:ngp
    for k=1:ngp
        f=quadpts(j); % Eta
        e=quadpts(k); % Zeta
        Xe=[X(1) Y(1) X(2) Y(2) X(3) Y(3) X(4) Y(4)];
        
        fi1=.25*(1-e)*(1-f);
        fi2=.25*(1+e)*(1-f);
        fi3=.25*(1+e)*(1+f);
        fi4=.25*(1-e)*(1+f);
        si1f =.25*[-(1-f) 0 (1-f) 0 (1+f) 0 -(1+f) 0]; %Si funtions
        si2f =.25*[0 -(1-f) 0 (1-f) 0 (1+f) 0 -(1+f)];
        si1e =.25*[-(1-e) 0 -(1+e) 0 (1+e) 0 (1-e) 0];
        si2e =.25*[0 -(1-e) 0 -(1+e) 0 (1+e) 0 (1-e)];
        J11=si1f*Xe';
        J12=si2f*Xe';
        J21=si1e*Xe';
        J22=si2e*Xe';
        J=[J11 J12; J21 J22]; %Jacobian
    A=1/det(J)*[J22 -J12 0 0; 0 0 -J21 J11;-J21 J11 J22 -J12]; 
    G=[si1f;si1e;si2f;si2e];
     r=fi1*X(1)+fi2*X(2)+fi3*X(3)+fi4*X(4);
    B1=A*G;
    B2=[fi1/r 0 fi2/r 0 fi3/r 0 fi4/r 0];
    B=[B1;B2];
   
    ke1=2*pi*B'*D*B*det(J)*weights(j)*weights(k)*r;
    ke(i).element=ke(i).element+ke1;
    
    %fb1=N'*det(J)*weights(j)*weights(k)*bodyforce; %Body Force
    %fb(i).element=fb(i).element+fb1;
   end
   end
        

  for j=1:max(size(D))*ndof
    for k=1:max(size(D))*ndof
       kg(conn(i,j),conn(i,k))=kg(conn(i,j),conn(i,k))+ke(i).element(j,k);
        
    end
     fg(conn(i,j),1)=fg(conn(i,j),1)+fet(i).element(j,1); %Traction force
    %fg(conn(i,j),1)=fg(conn(i,j),1)+fb(i).element(j,1); %Body Force
 end
 end



%% Applying Boundary Conditions and Solving
z=kg;
v=fg;
%Boundary Condition automatic
for i=1:max(size(bc))
z(bc(i),:)=0;z(:,bc(i))=0;z(bc(i),bc(i))=1;
end
u=linsolve(z,v); %displacement
r=kg*u-fg;       %reaction force  


%% Display
%fprintf('The element stiffness matrix for element %d is: \n',i);disp(ke((X)).element) %X - Element Number
%fprintf('The Global Stiffness matrix is: \n');disp(kg)  ;              %Global Stiffness Matrix
fprintf('The Displacement at point P is:\n');disp(u(9));disp(u(10));    %Dispacement at node
fprintf('The Reaction at point S is:\n');disp(r(1));disp(r(2));         %Reaction at node

%% Gauss Quadrature points generation
%Below function referenced from some blog to generate weights and absissas for
%any ngp(no of gauss points) given 
function[weights,quadpts]=gaussquads(ngp)
m=ngp-1;
syms x
P0=1;
P1=x;
for i=1:1:m
    Pn=((2.0*i+1)*x*P1-i*P0)/(i+1.0);
    P0=P1;
    P1=Pn;
end
if ngp==1
    Pn=P1;
end
Pn=expand(Pn);
quadpts=solve(vpa(Pn,32));
quadpts=sort(quadpts);
% Formula for weights is given at
% http://mathworld.wolfram.com/Legendre-GaussQuadrature.html Equation(13)
for k=1:1:ngp
    P0=1;
    P1=x;
    m=ngp;
    % Calculating P(n+1,x)
    for i=1:1:m
        Pn=((2.0*i+1)*x*P1-i*P0)/(i+1.0);
        P0=P1;
        P1=Pn;
    end
    Pn=P1;
    weights(k)=vpa(2*(1-quadpts(k)^2)/(ngp+1)^2/subs(Pn,x,quadpts(k))^2,32);
end
end




