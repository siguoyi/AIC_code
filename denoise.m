function [Pos_theta,theta_1s]=denoise(A,F,y)
theta_ls = (A(:,F)'*A(:,F))^(-1)*A(:,F)'*y;
        [val1,pos1]=sort(abs(theta_ls));%½µÐòÅÅÁÐ
        
       
        difference=diff(val1);

       difftol =(max(val1))/2^4;
       large = difference>difftol;

if any(large)
                pos = find(large,1,'first');
                
                p=pos+1;
else
                p=1;

end
        
       
        
%         c=1e-6;
%          for j=1:length(pos1)-1
%              if val1(j)-val1(j+1)>c
%                  c=val1(j)-val1(j+1);
%                  p=j;
%              end
%          end
        
        P = F(pos1(p:length(F)));
       Pos_theta=P;
       theta_1s = (A(:,P)'*A(:,P))^(-1)*A(:,P)'*y;