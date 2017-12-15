clear all;
global K M N long Fs J

%1epoch�̒���
long=1000;

%epoch�̌�
K=500;

%�M�����̐�
N=2;

%�T���v�����O���g��
Fs=8000;

%�ϑ��M���̐�
J=2;
plus=0;

%���̍����s��͈�l�����ŗ^����
% A=[1.0 0.5;0.5 1.0];
 
   %�܂��͈�l����
A=rand(N);

%�ǂݍ��މ���
[y,fs]=audioread('sc max fs8k CAN0001 130min.wav');
[y1,fs1]=audioread('sc max fs8k CAN1002 130min.wav');

% [ynew,fs2]=audioread('DIEMbunrishongoul.wav');

sum1=0;
sum2=0;
sum3=0;
sum4=0;
sum5=0;
sum6=0;

permutation=[0 1;1 0];

 tic

 for m=1:K
    
    for i=(m-1)*long+1:(m-1)*long+long
    
        for k=1:N
         
            %�ϑ��M���̌v�Z
  
            xshokichi(k,i) = A(k,1) * y(i,:) + A(k,2) * y1(i,:);
 
        end
          
        %�����M���̐���
        shuturyoku1(1,i)=xshokichi(1,i);
        shuturyoku2(1,i)=xshokichi(2,i);

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
filename = 'kongoushongoul.wav';
audiowrite(filename,shuturyoku1,Fs);

filename1 = 'kongoushongou2.wav';
audiowrite(filename1,shuturyoku2,Fs);


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
     
     %Step3�����s
     R1=reshape(max_eigvec,N,N);
     R2=reshape(max2_eigvec,N,N);
     
     
     M1hat=horzcat(R1(:,1),R2(:,1));
     M2hat=horzcat(R1(:,2),R2(:,2));
   
     Mlargehat=M1hat(:)*M1hat(:)'+M2hat(:)*M2hat(:)';
     
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
     
     Vnew=v2hat*inv(v1hat);
    
     [V,D] = eig(Vnew);
 
    G=inv(V)*A;
    
   if  trace(abs(permutation*G*G')) > trace(abs(G*G'))
         
       G = permutation*G;
       
     end
    
     gyousaidaichi=max(abs(G),[],2);
     retusaidaichi=max(abs(G),[],1);
     
     PIG=((G(1,1)^2+G(1,2)^2)/gyousaidaichi(1,1)^2-1)+((G(2,1)^2+G(2,2)^2)/gyousaidaichi(2,1)^2-1)+((G(1,1)^2+G(2,1)^2)/retusaidaichi(1,1)^2-1)+((G(1,2)^2+G(2,2)^2)/retusaidaichi(1,2)^2-1);
     PIGdb=10*log10(PIG);
     
    for m=1:K
    
         for i=(m-1)*long+1:(m-1)*long+long
     
                  shuturyoku3(1,i)=abs(G(1,1))*y(i,:) + abs(G(1,2))* y1(i,:);
                  shuturyoku4(1,i)=abs(G(2,1))*y(i,:) + abs(G(2,2))* y1(i,:);
    
            end
      
     end
    
     
     %wave�t�@�C����������    
filename2 = 'DIEMbunrishongoul.wav';
audiowrite(filename2,shuturyoku3,Fs);

filename3 = 'DIEMbunrishongou2.wav';
audiowrite(filename3,shuturyoku4,Fs);
    
      elapsedTime = toc
     
     
      for j=1:K*long
         
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
      Esir1=10*log10((sum1/K*long)/(sum2/K*long));
      
      Esir2=10*log10((sum4/K*long)/(sum3/K*long));
     
    
     heikin=(Esir1+Esir2)/2;

