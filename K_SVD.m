function [A,x]= K_SVD(y,codebook_size,errGoal)
%==============================
%input parameter
%   y - input signal
%   codebook_size - count of atoms
%output parameter
%   A - dictionary
%   x - coefficent
%reference:K-SVD:An Algorithm for Designing of Overcomplete Dictionaries
%          for Sparse Representation,Aharon M.,Elad M.etc
%==============================
if(size(y,2)<codebook_size)
    disp('codebook_size is too large or training samples is too small');
    return;
end
% initialization
[rows,cols]=size(y);
r=randperm(cols);
A=y(:,r(1:codebook_size));
A=A./repmat(sqrt(sum(A.^2,1)),rows,1);
ksvd_iter=10;
for k=1:ksvd_iter
    % sparse coding
    if nargin==2
        x=OMP(A,y,5.0/6*rows);
    elseif nargin==3
%         x=OMPerr(A,y,errGoal);
        x=[];
        for ii=1:1:rows
%           x=[x GoogleOMP(A,y(:,ii),errGoal)];
          x=[x GoogleOMP( y(:,ii),A,errGoal )];
        end
    end
    % update dictionary
    for m=1:codebook_size
        mindex=find(x(m,:));
        if ~isempty(mindex)
            mx=x(:,mindex);
            mx(m,:)=0;
            my=A*mx;
            resy=y(:,mindex);
            mE=resy-my;
            [u,s,v]=svds(mE,1);
            A(:,m)=u;
            x(m,mindex)=s*v';
        end
    end
end
