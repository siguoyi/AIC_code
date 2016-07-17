function [ hat_y ] = CS_OMP( y,A,K )

[A_rows,A_columns] = size(A);
    if A_rows<A_columns
        N=A_columns;
        M=A_rows;
    else
        N=A_rows;
        M=A_columns;
    end
hat_y=zeros(1,N);                                 %  待重构的谱域(变换域)向量                     
Aug_t=[];                                         %  增量矩阵(初始值为空矩阵)
r_n=y;                                            %  残差值

for times=1:K;                                    %  迭代次数(有噪声的情况下,该迭代次数为K)
    for col=1:N;                                  %  恢复矩阵的所有列向量
        product(col)=abs(A(:,col)'*r_n);          %  恢复矩阵的列向量和残差的投影系数(内积值) 
    end
    [val,pos]=max(product);                       %  最大投影系数对应的位置
    Aug_t=[Aug_t,A(:,pos)];                       %  矩阵扩充
    A(:,pos)=zeros(M,1);                          %  选中的列置零
    aug_y=(Aug_t'*Aug_t)^(-1)*Aug_t'*y;           %  最小二乘,使残差最小
    r_n=y-Aug_t*aug_y;                            %  残差
    pos_array(times)=pos;                         %  纪录最大投影系数的位置
end
hat_y(pos_array)=aug_y;                           %  重构的谱域向量
                    
end