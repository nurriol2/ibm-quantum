function [U,A,x,y,controls,num,V]=decomposition(n)
%n MUST BE GREATER THAN 1
%U  IS THE 2^n by 2^n UNITARY MATRIX TO BE DECOMPOSED
%A is the schemetable whose kth row indicate the type of control gate to be
%used in step k
%(x,y) is the vector of row and column indices (r(k),c(k)) is the kth entry
%to be annihilated
%num is the number of nontrivial gates used
%controls is the total number of controls
%V is a d by 4 array [V(i,1),V(i,2);V(i,3),V(i,4)] is the unitary matrix
%used in the ith gate

[x,y,A]=schemetable(n);

unitary= inputdlg('generate random unitary?[y/n]: ');
if strcmpi( unitary,'y') 
    U=randomunitary(n);
else 
    prompt = 'Please input 2^n by 2^n unitary matrix';
    U=input(prompt);
end
N=2^n;
d=N*(N-1)/2;
V=zeros(d,4);
Y=U;
controls=0;
num=d;
for j=1:d
[D,K,c]=ithgate(A(j,:),x(j,1),y(j,1),Y,n);
V(j,:)=K;
if V(j,:)==[1,0,0,1]
    num=num-1;
else
Y=D*Y;
controls=controls+c;
end
end
temp=0;
o=cellstr(strcat('(',num2str(x),',',num2str(y),')'))';
A=cellstr(A)';
for mm=1:N-1
    eval(['COLUMN', num2str(mm),'=','[o(temp+1:temp+2^n-mm);A(temp+1:temp+2^n-mm)]'])
    for kk=1:2^n-mm
        eval(['Vr',num2str(x(temp+kk)),'c',num2str(y(temp+kk)),'=','[V(temp+kk,1:2);V(temp+kk,3:4)]'])
    end
    temp=temp+N-mm;
end

function [W]=randomunitary(n) 
%THIS FUNCTION GENERATES A RANDOM 2^n by 2^n UNITARY MATRIX
W=rand(2^n);
H=0.5*(W+W');
W=expm(1i*H);

function [x,y,A] = schemetable(n) 
%THIS FUNCTION GENERATES THE SCHEME TABLE FOR n

N=2^(n);              %dimension of matrix
d=N*(N-1)/2;          %number of gates used 

%ORDER OF ANNIHILATION, COLUMN INDEX
y=zeros(d,1);         
temp1=0;
  for j=1:N-1;
      y(temp1+1:temp1+N-j,1)=j;
      temp1=temp1+2^n-j;
  end


%ORDER OF ANNIHILATION, ROW INDEX AND GATES

%initializations
x=zeros(d,1);          
A=repmat('*',d,n);
x(1,1)=2;         
A(1,n)='T';       


for k=2:n           %loop index signifies leading 2^k by 2^k subblock
    temp2=2^(k-1);  
    
    %first column entries and gates 
    x(temp2:2*temp2-1,1)=[x(1:temp2-1,1)+temp2; temp2+1]; %column 1, lower half of 2^k leading column
    A(temp2:2*temp2-1,n-k+1:n)=[A(1:temp2-1,n-k+1:n); ['T',repmat('*',1,k-1);]];
    for i=1:k-1
        A(temp2+2^(i)-2,n-k+1)='1'; %copies control free gate that annihilates 2^i+1 and adds a 1-control
    end
    
    temp3=2^n-1;   %this temp variable becomes the sum of 2^n-(j-1), where j is the previous column dealt with
    for j=2:temp2  %loop adapts first column, lower left of 2^k block to columns indexed by j
        x(temp3+temp2-j+1:temp3+2*temp2-j,1)=Fell(k,j,x(temp2:2*temp2-1,1));
        A(temp3+temp2-j+1:temp3+2*temp2-j,n-k+1:n)=Gell(A(temp2:2*temp2-1,n-k+1:n),j);
       temp3=temp3+N-j;
    end
    
    
    for jj=1:temp2-1 %loop upper left [and adds 2^(k-1)] to lower right of 2^k block
        bb=(jj-1)*(N-jj/2);
        x(temp3+1:temp3+temp2-jj,1)=x(bb+1:bb+temp2-jj,1)+temp2;
        A(temp3+1:temp3+temp2-jj,n-k+1:n)=[repmat('1',temp2-jj,1),A(bb+1:bb+temp2-jj,n-k+2:n)];
        temp3=temp3+N-temp2-jj;
    end
end

function [D,K,c]=ithgate(Ai,xi,yi,U,n)
c=0;
D=1;
if yi==2^n-1
    D=U';
    K=[conj(U(2^n-1,2^n-1));conj(U(2^n,2^n-1));conj(U(2^n-1,2^n));conj(U(2^n,2^n))];
    c=c+n-1;
else
    for k=n:-1:1
        if isequal(Ai(k),'0')==1
            D=kron([1,0;0,0],D);
            c=c+1;
        elseif isequal(Ai(k),'1')==1
            D=kron([0,0;0,1],D);
            c=c+1;
        elseif isequal(Ai(k),'T')==1 && U(xi,yi)==0
             K=[1,0,0,1];
             D=kron(zeros(2,2),D);
        elseif isequal(Ai(k),'T')==1 && (bitget(yi-1,n-k+1)==0)
            a=U(Fell(n,2^(n-k)+1,xi),yi);
            b=U(xi,yi);
            z=sqrt(abs(a)^2+abs(b)^2);
            K=(1/z)*[a,-b,conj(b),conj(a)];
            D=kron(-1*eye(2,2)+[K(1,1:2);K(1,3:4)],D);
        elseif isequal(Ai(k),'T')==1 && (bitget(yi-1,n-k+1)==1)
            a=U(Fell(n,2^(n-k)+1,xi),yi);
            b=U(xi,yi);
            z=sqrt(abs(a)^2+abs(b)^2);
            K=(1/z)*[conj(a),conj(b),-b,a];
            D=kron(-1*eye(2,2)+[K(1,1:2);K(1,3:4)],D);
        else
            D=kron(eye(2,2),D);
        end
    end
    D=D+eye(2^n,2^n);
end

function [Y]=Gell(X,l)
%l is an integer from 1 to 2^k-1
%X must be p by k
%Y=G_l(X)
Y=X;
[p,k]=size(X);
C=repmat('1',1,k);
s=dec2bin(l-1);
r=size(s,2);
Y(p,1)='T';
  for m=1:r
    if bitget(l-1,m)==1
        for t=1:p-1
            if X(t,k-m+1)=='1'
            Y(t,k-m+1)='0';
            end  
        end
    Y(p,k-m+1)='1';
    end
  end
  for t=1:p-1
     if size(intersect(X(t,1:k-r),C),2)==0
     Y(t,1)='1';
     end 
  end

function [ v ] = Fell(n,r,u)
%Fell takes a vector of integers u and sends it to the vector of integers u, the binary
%representation of u(i,j) and v(i,j) differ precisely in places where the binary
%digit of r (in a word of length n) is 1 
ub=dec2bin(u(:,1)-1,n);
rv=r*ones(size(u));
rb=dec2bin(rv(:,1)-1,n);
flip=mod(ub+rb,2);
fbits=cellstr(num2str(flip));
v=bin2dec(fbits(:,1))+1;


