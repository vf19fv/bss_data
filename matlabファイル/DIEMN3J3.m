clear all;
global K M N long Fs J

%1epochの長さ
long=10000;

%epochの個数
K=1000;

%信号源の数
N=3;

%観測信号の数
J=3;

%サンプリング周波数
Fs=8000;

plus=0;

%元の混合行列は一様乱数で与える
A=[1.0 0.5 0.25;0.5 1.0 0.5;0.25 0.5 1.0];

%混合行列Aを乱数
%A=rand(K);
    
%読み込む音声
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
         
         %観測信号の計算
  
         xshokichi(k,i) = A(k,1) * y(i,:) + A(k,2) * y1(i,:)+A(k,3)*y2ongen(i,:);
 
          end
          
       %混合信号の生成
       shuturyoku1(1,i)=xshokichi(1,i);
       shuturyoku2(1,i)=xshokichi(2,i);
       shuturyoku3(1,i)=xshokichi(3,i);

       Tk{i} = xshokichi(:,i)*xshokichi(:,i)';
     
     Tn{1} = 0; 
     Tn{i+1} = Tn{i} + Tk{i};
     
   end
   
     %quasistationary性を仮定した共分散行列を計算
     Tnew{1,m} = Tn{(m-1)*long+long+1}; 
     
     %ベクタライズ
     yvec{1,m}=Tnew{1,m}(:);
     
     %共分散行列を計算
     yseki{1,m}=yvec{1,m}*yvec{1,m}';
     
     %和を計算
     plus=plus+yseki{1,m};
 
     
 end

       %waveファイル書き込み    
     filename = 'kongoushongouN3J3-l.wav';
     audiowrite(filename,shuturyoku1,Fs);

     filename1 = 'kongoushongouN3J3-2.wav';
     audiowrite(filename1,shuturyoku2,Fs);
     
     filename2 = 'kongoushongouN3J3-3.wav';
     audiowrite(filename2,shuturyoku3,Fs); 
     
     
%固有値分解
[Vfai,Dfai] = eig(plus);

%N^2個の固有値からK個の大きい固有値に対応する固有ベクトルを取り出す
 for i=1:J^2
         
         Dfainew(1,i)=Dfai(i,i);
         
 end
      
     %最大となるSの要素が何番目にあるのか探索
     [xfainew,yfainew] = max(Dfainew);
    
     %最大固有値に対応する固有ベクトルを求める
     max_eigvec = Vfai(:, yfainew);
    
     %最大固有値,それに対応する固有ベクトルを削除
     Dfai(:, yfainew) = [];
     Vfai(:, yfainew) = [];
     
     
     [x2fainew,y2fainew] = max(max(Dfai));
     max2_eigvec = Vfai(:, y2fainew);
     
     
     %3番目に大きい固有値,それに対応する固有ベクトルを削除
     Dfai(:, y2fainew) = [];
     Vfai(:, y2fainew) = [];
     
     [x3fainew,y3fainew] = max(max(Dfai));
     max3_eigvec = Vfai(:, y3fainew);
     
     %Step3を実行
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
  
%      %最大となるSの要素が何番目にあるのか探索
      [xnew,ynew] = max(D6new);
       
     %Step6において最大固有値に対応する固有ベクトルを求める
      max6_eigvec = V6(:, ynew);
    
     %最大固有値,それに対応する固有ベクトルを削除
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
     
     %waveファイル書き込み    
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

