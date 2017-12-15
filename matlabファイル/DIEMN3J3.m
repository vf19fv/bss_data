clear all;
global K M N long Fs J

%1epoch�̒���
long=10000;

%epoch�̌�
K=1000;

%�M�����̐�
N=3;

%�ϑ��M���̐�
J=3;

%�T���v�����O���g��
Fs=8000;

plus=0;

%���̍����s��͈�l�����ŗ^����
A=[1.0 0.5 0.25;0.5 1.0 0.5;0.25 0.5 1.0];

%�����s��A�𗐐�
%A=rand(K);
    
%�ǂݍ��މ���
[y,fs]=audioread('sc max fs8k CAN0001 130min.wav');
[y1,fs1]=audioread('sc max fs8k CAN1002 130min.wav');
[y2ongen,fs2]=audioread('sc max fs8k CAN1001 130min.wav');

permutation1=[0 1 0;1 0 0;0 0 1];

permutation2=[1 0 0;0 0 1;0 1 0];

permutation3=[0 0 1;0 1 0;1 0 0];

permutation4=[0 0 1;1 0 0;0 1 0];

permutation5=[0 1 0;0 0 1;1 0 0];

sum1=0;
sum2=0;
sum3=0;
sum4=0;
sum5=0;
sum6=0;


 tic

 for m=1:K
    
      for i=(m-1)*long+1:(m-1)*long+long
    
          for k=1:J
         
         %�ϑ��M���̌v�Z
  
         xshokichi(k,i) = A(k,1) * y(i,:) + A(k,2) * y1(i,:)+A(k,3)*y2ongen(i,:);
 
          end
          
       %�����M���̐���
       shuturyoku1(1,i)=xshokichi(1,i);
       shuturyoku2(1,i)=xshokichi(2,i);
       shuturyoku3(1,i)=xshokichi(3,i);

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
 
     
 end

       %wave�t�@�C����������    
     filename = 'kongoushongouN3J3-l.wav';
     audiowrite(filename,shuturyoku1,Fs);

     filename1 = 'kongoushongouN3J3-2.wav';
     audiowrite(filename1,shuturyoku2,Fs);
     
     filename2 = 'kongoushongouN3J3-3.wav';
     audiowrite(filename2,shuturyoku3,Fs); 
     
     
%�ŗL�l����
[Vfai,Dfai] = eig(plus);

%N^2�̌ŗL�l����K�̑傫���ŗL�l�ɑΉ�����ŗL�x�N�g�������o��
 for i=1:J^2
         
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
     
     
     %3�Ԗڂɑ傫���ŗL�l,����ɑΉ�����ŗL�x�N�g�����폜
     Dfai(:, y2fainew) = [];
     Vfai(:, y2fainew) = [];
     
     [x3fainew,y3fainew] = max(max(Dfai));
     max3_eigvec = Vfai(:, y3fainew);
     
     %Step3�����s
     R1=reshape(max_eigvec,J,J);
     R2=reshape(max2_eigvec,J,J);
     R3=reshape(max3_eigvec,J,J);
     
     M1hat=horzcat(R1(:,1),R2(:,1),R3(:,1));
     M2hat=horzcat(R1(:,2),R2(:,2),R3(:,2));
     M3hat=horzcat(R1(:,3),R2(:,3),R3(:,3));
   
     Mlargehat=M1hat(:)*M1hat(:)'+M2hat(:)*M2hat(:)'+M3hat(:)*M3hat(:)';
     
     [V6,D6] = eig(Mlargehat);
     
      for i=1:J^2
          
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
    
     v1hat=reshape(max6_eigvec,J,J);
     v2hat=reshape(max62_eigvec,J,J);
   
    
     [Vi,Di] = eig(v2hat*inv(v1hat));
 
     G=inv(Vi)*A;
     
     if  trace(abs(permutation1*G*G')) > trace(abs(G*G'))
         
       G = permutation1*G;
       
     else if trace(abs(permutation2*G*G')) > trace(abs(G*G'))
             
       G = permutation2*G;
       
     else if trace(abs(permutation3*G*G')) > trace(abs(G*G'))
             
       G = permutation3*G;
       
     else if trace(abs(permutation4*G*G')) > trace(abs(G*G'))
             
       G = permutation4*G;
       
     else if trace(abs(permutation5*G*G')) > trace(abs(G*G'))
             
       G = permutation5*G;  
       
         end
     
         end
         end
         
         end
         
     end
    
     gyousaidaichi=max(abs(G),[],2);
     retusaidaichi=max(abs(G),[],1);
    
     PIGgyou=((G(1,1)^2+G(1,2)^2+G(1,3)^2)/gyousaidaichi(1,1)^2-1)+((G(2,1)^2+G(2,2)^2+G(2,3)^2)/gyousaidaichi(2,1)^2-1)+((G(3,1)^2+G(3,2)^2+G(3,3)^2)/gyousaidaichi(3,1)^2-1);

     PIGretu=((G(1,1)^2+G(2,1)^2+G(3,1)^2)/retusaidaichi(1,1)^2-1)+((G(1,2)^2+G(2,2)^2+G(3,2)^2)/retusaidaichi(1,2)^2-1)+((G(1,3)^2+G(2,3)^2+G(3,3)^2)/retusaidaichi(1,3)^2-1);
      
      po=(PIGgyou+PIGretu)/6;
      
    PIGdb=10*log10(po);
   
    
        for m=1:K
    
         for i=(m-1)*long+1:(m-1)*long+long
     
                  shuturyoku3(1,i)=abs(G(1,1))*y(i,:) + abs(G(1,2))* y1(i,:)+abs(G(1,3))*y2ongen(i,:);
                  shuturyoku4(1,i)=abs(G(2,1))*y(i,:) + abs(G(2,2))* y1(i,:)++abs(G(2,3))*y2ongen(i,:);
                  shuturyoku5(1,i)=abs(G(3,1))*y(i,:) + abs(G(3,2))* y1(i,:)++abs(G(3,3))*y2ongen(i,:);
    
            end
      
     end
     
     %wave�t�@�C����������    
filename1 = 'DIEMbunrishongouN3J3-l.wav';
audiowrite(filename1,shuturyoku3,Fs);

filename2 = 'DIEMbunrishongouN3J3-2.wav';
audiowrite(filename2,shuturyoku4,Fs);

filename3 = 'DIEMbunrishongouN3J3-3.wav';
audiowrite(filename3,shuturyoku5,Fs);

      elapsedTime = toc
     
      
  for j=1:K*long
         
         nijyou(j,1)=(G(1,1)*y(j,:)).^2;
         nijyou2(j,1)=(G(1,2)*y1(j,:)).^2;
         nijyou3(j,1)=(G(1,3)*y2ongen(j,:)).^2;
          
           sum1=sum1+nijyou(j,1);
           sum2=sum2+nijyou2(j,1)+nijyou3(j,1);
          
     
% %           sum1=sum1+nijyou;
%           
           nijyou4(j,1)=(G(2,2)*y1(j,:)).^2;
           nijyou5(j,1)=(G(2,1)*y(j,:)).^2;
           nijyou6(j,1)=(G(2,3)*y2ongen(j,:)).^2;
           
           sum3=sum3+nijyou4(j,1);
           sum4=sum4+nijyou5(j,1)+nijyou6(j,1);
           
           nijyou7(j,1)=(G(3,3)*y2ongen(j,:)).^2;
           nijyou8(j,1)=(G(3,1)*y(j,:)).^2;
           nijyou9(j,1)=(G(3,2)*y1(j,:)).^2;
           
           sum5=sum5+nijyou7(j,1);
           sum6=sum6+nijyou8(j,1)+nijyou9(j,1);
          
%           sum3=sum3+nijyou3(j,1);
%           sum4=sum4+nijyou4(j,1);
%      
       end
          
        Esir1=10*log10((sum1/K*long)/(sum2/K*long));
      
        Esir2=10*log10((sum3/K*long)/(sum4/K*long));
%      
       Esir3=10*log10((sum5/K*long)/(sum6/K*long));
%     
       heikin=(Esir1+Esir2+Esir3)/K;

