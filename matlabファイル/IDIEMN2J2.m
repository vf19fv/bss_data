clear all;
global K M N long

long=10000;
M=100;
K=2;
N=2;
plus=0;

%���̍����s��͈�l�����ŗ^����
    A=[1.0 0.5;0.5 1.0];
%   A=rand(K);

[y,fs]=audioread('sc max fs8k CAN0001 130min.wav');
[y1,fs1]=audioread('sc max fs8k CAN1002 130min.wav');

sum1=0;
sum2=0;
sum3=0;
sum4=0;
sum5=0;
sum6=0;

 tic

 for m=1:M   
    
      for i=(m-1)*long+1:(m-1)*long+long
    
          for k=1:K
         
         %�ϑ��M���̌v�Z
  
         xshokichi(k,i) = A(k,1) * y(i,:) + A(k,2) * y1(i,:);
 
       end
          
%           xnew(1,i)=xnoisefree(1,i)+noise1(i,1);
%           xnew(2,i)=xnoisefree(2,i)+noise2(i,1);

       Tk{i} = xshokichi(:,i)*xshokichi(:,i)';
     
     Tn{1} = 0; 
     Tn{i+1} = Tn{i} + Tk{i};
     
   end
   
     %quasistationary�������肵�������U�s����v�Z
     Tnew{1,m} = Tn{(m-1)*long+long+1}; 
     
     %�x�N�^���C�Y
     yvec{1,m}=Tnew{1,m}(:);
     
     %�����U�s����v�Z
     yseki{1,m}=yvec{1,m}*yvec{1,m}';
     
     %�a���v�Z
     plus=plus+yseki{1,m};
   
    % z{1,m}=Tnew{1,m}(:)./trace(Tnew{1,m});
     
     %�S�ăx�N�^���C�Y
    % znew{1,m}=z{1,m}(:);
     
     %�S�Ă̗v�f��2��
    % znormu{1,m}=power(znew{1,m},2);
     
     %2�m�����̒l���v�Z
    % S(1,m)=sum(znormu{1,m});
     
end


%�ŗL�l����
[Vfai,Dfai] = eig(plus);

%N^2�̌ŗL�l����K�̑傫���ŗL�l�ɑΉ�����ŗL�x�N�g�������o��
 for i=1:K^2
         
         Dfainew(1,i)=Dfai(i,i);
         
 end
      
     %�ő�ƂȂ�S�̗v�f�����Ԗڂɂ���̂��T��
     [xfainew,yfainew] = max(Dfainew);
    
     %�ő�ŗL�l�ɑΉ�����ŗL�x�N�g�������߂�
     max_eigvec = Vfai(:, yfainew);
    
     %�ő�ŗL�l,����ɑΉ�����ŗL�x�N�g�����폜
     Dfai(:, yfainew) = [];
     Vfai(:, yfainew) = [];
     
     [x2fainew,y2fainew] = max(max(Dfai));
     max2_eigvec = Vfai(:, y2fainew);
     
     %Step3�����s
     R1=reshape(max_eigvec,K,K);
     R2=reshape(max2_eigvec,K,K);
     
     
     M1hat=horzcat(R1(:,1),R2(:,1));
     M2hat=horzcat(R1(:,2),R2(:,2));
   
     %�����܂ł�DIEM�Ɠ���
     Mlargehat=M1hat(:)*M1hat(:)'+M2hat(:)*M2hat(:)';
     
     [V6,D6] = eig(Mlargehat);
     
      for i=1:K^2
          
        D6new(1,i)=D6(i,i);
          
      end
  
%      %�ő�ƂȂ�S�̗v�f�����Ԗڂɂ���̂��T��
      [xnew,ynew] = max(D6new);
       
     %Step6�ɂ����čő�ŗL�l�ɑΉ�����ŗL�x�N�g�������߂�
      max6_eigvec = V6(:, ynew);
    
     %�ő�ŗL�l,����ɑΉ�����ŗL�x�N�g�����폜
    D6(:, ynew) = [];
    V6(:, ynew) = [];
    
    [x62fainew,y62fainew] = max(max(D6));
    
    max62_eigvec = V6(:, y62fainew);
    
    
     v1hat=reshape(max6_eigvec,K,K);
     v2hat=reshape(max62_eigvec,K,K);
     
     Vnew=inv(v1hat)*v2hat;
     
     
     
    %�ŗL�l�s��V
     [E,D] = eig(Vnew);
 
     lambda1=v1hat*E(:,1);
     lambda2=v2hat*E(:,1);
      Fai=horzcat(lambda1,lambda2);
      
      maxlambda=Fai*Fai';
      
      
      [Vki,D1] = eig(maxlambda);

%N^2�̌ŗL�l����K�̑傫���ŗL�l�ɑΉ�����ŗL�x�N�g�������o��
 for i=1:K         
     
         Dfai1(1,i)=D1(i,i);
         
 end
      
     %�ő�ƂȂ�S�̗v�f�����Ԗڂɂ���̂��T��
     [xfai1,yfai1] = max(Dfai1);
    
     %�ő�ŗL�l�ɑΉ�����ŗL�x�N�g�������߂�
     max_eigvecpart1 = Vki(:, yfai1);
      
     
  
     lambda3=v1hat*E(:,2);
     lambda4=v2hat*E(:,2);
     
     Fai2=horzcat(lambda3,lambda4);
     
      maxlambda2=Fai2*Fai2';
      
      [V2k,D2k] = eig(maxlambda2);
      
       for i=1:K      
           
         Dfai2(1,i)=D2k(i,i);
         
       end
       
     %�ő�ƂȂ�S�̗v�f�����Ԗڂɂ���̂��T��
     [xfai2,yfai2] = max(Dfai2);
 
     %�ő�ŗL�l�ɑΉ�����ŗL�x�N�g�������߂�
     max_eigvecpart2 = V2k(:, yfai2);
     
     
     
     Vsa=horzcat(max_eigvecpart1,max_eigvecpart2);
      
    G=inv(Vsa)*A;
    
    
     gyousaidaichi=max(abs(G),[],2);
     retusaidaichi=max(abs(G),[],1);
%     
     
      PIG=((G(1,1)^2+G(1,2)^2)/gyousaidaichi(1,1)^2-1)+((G(2,1)^2+G(2,2)^2)/gyousaidaichi(2,1)^2-1)+((G(1,1)^2+G(2,1)^2)/retusaidaichi(1,1)^2-1)+((G(1,2)^2+G(2,2)^2)/retusaidaichi(1,2)^2-1);
     PIGdb=10*log10(PIG);
     
            
      elapsedTime = toc
     
     
      for j=1:M*long
         
           nijyou(j,1)=(G(1,1)*y(j,:))^2;
           nijyou2(j,1)=(G(1,2)*y1(j,:))^2;
         
            sum1=sum1+nijyou(j,1);
            sum2=sum2+nijyou2(j,1);
%          
%            sum1=sum1+nijyou;
          
           nijyou3(j,1)=(G(2,1)*y(j,:))^2;
           nijyou4(j,1)=(G(2,2)*y1(j,:))^2;
         
            sum3=sum3+nijyou3(j,1);
            sum4=sum4+nijyou4(j,1);
     
      end
      
      
       
      %�Ίp�v�f����Ίp�v�f
      Esir1=10*log10((sum1/M*long)/(sum2/M*long));
      
      Esir2=10*log10((sum4/M*long)/(sum3/M*long));
     
    
     heikin=(Esir1+Esir2)/2;

