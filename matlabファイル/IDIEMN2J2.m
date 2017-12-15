clear all;
global K M N long

long=10000;
M=100;
K=2;
N=2;
plus=0;

%元の混合行列は一様乱数で与える
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
         
         %観測信号の計算
  
         xshokichi(k,i) = A(k,1) * y(i,:) + A(k,2) * y1(i,:);
 
       end
          
%           xnew(1,i)=xnoisefree(1,i)+noise1(i,1);
%           xnew(2,i)=xnoisefree(2,i)+noise2(i,1);

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
   
    % z{1,m}=Tnew{1,m}(:)./trace(Tnew{1,m});
     
     %全てベクタライズ
    % znew{1,m}=z{1,m}(:);
     
     %全ての要素を2乗
    % znormu{1,m}=power(znew{1,m},2);
     
     %2ノルムの値を計算
    % S(1,m)=sum(znormu{1,m});
     
end


%固有値分解
[Vfai,Dfai] = eig(plus);

%N^2個の固有値からK個の大きい固有値に対応する固有ベクトルを取り出す
 for i=1:K^2
         
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
     
     %Step3を実行
     R1=reshape(max_eigvec,K,K);
     R2=reshape(max2_eigvec,K,K);
     
     
     M1hat=horzcat(R1(:,1),R2(:,1));
     M2hat=horzcat(R1(:,2),R2(:,2));
   
     %ここまではDIEMと同じ
     Mlargehat=M1hat(:)*M1hat(:)'+M2hat(:)*M2hat(:)';
     
     [V6,D6] = eig(Mlargehat);
     
      for i=1:K^2
          
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
    
    
     v1hat=reshape(max6_eigvec,K,K);
     v2hat=reshape(max62_eigvec,K,K);
     
     Vnew=inv(v1hat)*v2hat;
     
     
     
    %固有値行列がV
     [E,D] = eig(Vnew);
 
     lambda1=v1hat*E(:,1);
     lambda2=v2hat*E(:,1);
      Fai=horzcat(lambda1,lambda2);
      
      maxlambda=Fai*Fai';
      
      
      [Vki,D1] = eig(maxlambda);

%N^2個の固有値からK個の大きい固有値に対応する固有ベクトルを取り出す
 for i=1:K         
     
         Dfai1(1,i)=D1(i,i);
         
 end
      
     %最大となるSの要素が何番目にあるのか探索
     [xfai1,yfai1] = max(Dfai1);
    
     %最大固有値に対応する固有ベクトルを求める
     max_eigvecpart1 = Vki(:, yfai1);
      
     
  
     lambda3=v1hat*E(:,2);
     lambda4=v2hat*E(:,2);
     
     Fai2=horzcat(lambda3,lambda4);
     
      maxlambda2=Fai2*Fai2';
      
      [V2k,D2k] = eig(maxlambda2);
      
       for i=1:K      
           
         Dfai2(1,i)=D2k(i,i);
         
       end
       
     %最大となるSの要素が何番目にあるのか探索
     [xfai2,yfai2] = max(Dfai2);
 
     %最大固有値に対応する固有ベクトルを求める
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
      
      
       
      %対角要素÷非対角要素
      Esir1=10*log10((sum1/M*long)/(sum2/M*long));
      
      Esir2=10*log10((sum4/M*long)/(sum3/M*long));
     
    
     heikin=(Esir1+Esir2)/2;

