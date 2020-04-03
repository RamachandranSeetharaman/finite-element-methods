clc
clear all
 
%% Data
E=200e9;            %Youngs Modulus
mew=0.3;            %Poisson Ratio
h=.02;              %thickness
ndof=2;             %degree of freedom per node
endcoo=[0 0;4 0;4 4;0 4 ]; %End coordinates of each node
nnx=10; % No of nodes in X direction
nny=5; % No of nodes in Y direction




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
ne=2*2*(nex+ney);       % Total No of elements
fg=zeros(1,tdof);       %force matrix initialization
kg=zeros(tdof,tdof);    %Global Elemental stiffness matrix generation

D=(E/(1-mew^2))*[1 mew 0;mew 1 0; 0 0 (1-mew)/2];

%% Force Input and Boundary Conditions
%Force inout
fg(49)=10000;fg(39)=20000;fg(29)=20000;fg(19)=20000;fg(9)=10000; %Force input in nodes
%Input force corresponding to coordinate system. 
%eg., Downward force is taken negative
%bodyforce=[100;200]; %bodyforce [bx;by];

%Boundary Condition
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
 for i=1:size(quads);
     for j=1:nnpe
         X(j)=Qx(quads(i,j));
         Y(j)=Qy(quads(i,j));
     end
     x1=X(1);x2=X(2);x3=X(3);x4=X(4);y1=Y(1);y2=Y(2);y3=Y(3);y4=Y(4);
    xlim([0 2*lx])
    ylim([0 2*ly])
    plot(X,Y)
    hold on
     %Ae=polyarea(X,Y)
    A=[1 x1 y1;1 x2 y2;1 x3 y3];
    Ae=.5*det(A);
    x12=x1-x2;x23=x2-x3;x13=x1-x3;x31=-x13;x32=-x23;x21=-x12;
    y12=y1-y2;y23=y2-y3;y13=y1-y3;y31=-y13;y32=-y23;y21=-y12;
    B(i).element=(1/(2*Ae))*[y23 0 y31 0 y12 0;0 x32 0 x13 0 x21;x32 y23 x13 y31 x21 y12];
    ke(i).element=B(i).element'*D*B(i).element*Ae*h;%Elemental Stiffness Matrix
       %fb1=N'*det(J)*weights(j)*weights(k)*bodyforce; %Body Force
    %fb(i).element=fb(i).element+fb1;
   
  
        
  %fprintf('The element stiffness matrix for element %d is: \n',i);disp(ke(i).element)
  for j=1:max(size(D))*ndof %8 since node per element is 4 and each node has 2 dof
    for k=1:(size(D))*ndof
       kg(conn(i,j),conn(i,k))=kg(conn(i,j),conn(i,k))+ke(i).element(j,k);
        
    end
       %fg(conn(i,j),1)=fg(conn(i,j),1)+fb(i).element(j,1); %Body Force
  end
end
hold off

    

%% Applying Boundary Conditions and Solving
z=kg;
v=fg;
%Boundary Condition automatic
for i=1:max(size(bc))
z(bc(i),:)=0;z(:,bc(i))=0;z(bc(i),bc(i))=1;
end
u=linsolve(z,v'); %displacement
r=kg*u-fg';       %reaction force  


%% Display
%fprintf('The element stiffness matrix for element %d is: \n',i);disp(ke((X)).element) %X - Element Number
%fprintf('The Global Stiffness matrix is: \n');disp(kg)  ;              %Global Stiffness Matrix
fprintf('The Displacement at point P is:\n');disp(u(9));disp(u(10));    %Dispacement at node
fprintf('The Reaction at point S is:\n');disp(r(1));disp(r(2));         %Reaction at node






