%% 適応制御　6章の例題
% 例6.1 バックステッピング法
clear
% close all
%% システム定義
s=tf('s');
tf_y=1/(s^2+1.4*s+1); % 制御対象のモデル(実際には我々は係数を知らない)
tf_ym=1/(s+1)^2; % 規範モデル
Ws=(s+1)^2; % 規範モデルの分母
nn=2; % システムの相対次数（制御対象の相対次数は既知としている）
n=2; % システムの次数（制御対象の次数も既知としている）

%% 制御対象の実現
ss_y=ss(tf_y);
ss_Ws=ss(1/Ws);

%% 状態変数フィルタ（設計変数）
% (F,g)が可制御になるように選ぶ
F=[-1 0;
    0 -2];
g=[1;1];

%% その他必要な定数
b0=ss_y.c*ss_y.A^(nn-1)*ss_y.b;
lambda=1; % >0

%% 適応ゲイン 今回はn次元のフィルタを用いた
G=10*eye(5);

%% 入力
dt=0.025; % 刻み時間はマシンスペックとか計算誤差を見ながら決めてます。
t=0:dt:40;
r=sin(t);

%% 計算結果を入れる箱
z1=zeros(1,length(t)); % 出力誤差
z2=zeros(1,length(t)); % uf1とalpha1の誤差
ym=zeros(2,length(t)); % 規範モデルの応答
ym_true=zeros(1,length(t)); % 規範モデルの応答
dym=zeros(2,length(t)); % 規範モデルの応答の微分
dym_true=zeros(1,length(t)); % 規範モデルの応答の微分
ddym_true=zeros(1,length(t)); % 規範モデルの応答の微分
dYm=zeros(1,length(t)); % 規範モデルの応答のとその微分の和（計算で使うんです）
Ym=zeros(1,length(t)); % 規範モデルの応答のとその微分の和（計算で使うんです）
y=zeros(2,length(t)); % 実際の制御対象の応答
y_true=zeros(1,length(t)); % 実際の制御対象の応答
yf1=zeros(1,length(t)); % 実際の制御対象の応答を積分して得る値
dyf1=zeros(1,length(t)); % 実際の制御対象の応答を積分して得る値
u=zeros(1,length(t)); % 制御対象への入力
uf1=zeros(1,length(t)); % 入力を積分して得る信号
theta=zeros(2*n+1,length(t)); % 適応
omega1=zeros(2*n+1,length(t)); % 適応や入力に使う変数群
domega1=zeros(2*n+1,length(t)); % 適応や入力に使う変数群
ep=zeros(2,length(t)); % 指数減衰項らしいけどよくわからない
ep1=zeros(1,length(t)); % 指数減衰項の第一成分
alpha1=zeros(1,length(t)); % 入力を積分して得る信号を再現するための信号
p=zeros(1,length(t)); % alphaを生成するために必要な同定変数
fai=zeros(1,length(t)); % alphaを生成するために必要な同定変数
v1=zeros(2,length(t)); % 状態変数フィルタの状態量
v2=zeros(2,length(t)); % 状態変数フィルタの状態量
dv1=zeros(2,length(t)); % 状態変数フィルタの状態量
dv2=zeros(2,length(t)); % 状態変数フィルタの状態量
b0_hat=zeros(1,length(t));
tauth1=zeros(5,length(t));
tauth2=zeros(5,length(t));
taub2=zeros(1,length(t));
alpha2=zeros(1,length(t));
ganma1=zeros(1,length(t));
ganmath1=zeros(5,length(t));
beta2=zeros(1,length(t));
%% 計算
c1=3;
c2=3;
d1=3;
d2=3;
g1=5;
g3=5;
tic
for ii=1:length(t)
   %% ステップ1
    y_true(ii)=ss_y.c*y(:,ii); % yは状態量で，実際の出力をy_trueに設置
    % yf1を求める
    K1=dt*(-lambda*yf1(ii)+y_true(ii));
    K2=dt*(-lambda*(yf1(ii)+0.5*K1)+y_true(ii));
    K3=dt*(-lambda*(yf1(ii)+0.5*K2)+y_true(ii));
    K4=dt*(-lambda*(yf1(ii)+K3)+y_true(ii));
    yf1(ii+1)=yf1(ii)+(K1+2*K2+2*K3+K4)/6;
    dyf1(:,ii)=(yf1(:,ii+1)-yf1(:,ii))/dt;
    
    % v1(ii+1)を求める
    K1=dt*(g*yf1(ii)+F*v1(:,ii));
    K2=dt*(g*yf1(ii)+F*(v1(:,ii)+0.5*K1));
    K3=dt*(g*yf1(ii)+F*(v1(:,ii)+0.5*K2));
    K4=dt*(g*yf1(ii)+F*(v1(:,ii)+K3));
    v1(:,ii+1)=v1(:,ii)+(K1+2*K2+2*K3+K4)/6;
    dv1(:,ii)=(v1(:,ii+1)-v1(:,ii))/dt;
    
    % v2(ii+1)を求める
    K1=dt*(g*uf1(ii)+F*v2(:,ii));
    K2=dt*(g*uf1(ii)+F*(v2(:,ii)+0.5*K1));
    K3=dt*(g*uf1(ii)+F*(v2(:,ii)+0.5*K2));
    K4=dt*(g*uf1(ii)+F*(v2(:,ii)+K3));
    v2(:,ii+1)=v2(:,ii)+(K1+2*K2+2*K3+K4)/6;
    dv2(:,ii)=(v2(:,ii+1)-v2(:,ii))/dt;
    
   % ym(ii+1)を求める
    K1=dt*(ss_Ws.a*ym(:,ii)+ss_Ws.b*r(ii));
    K2=dt*(ss_Ws.a*(ym(:,ii)+0.5*K1)+ss_Ws.b*r(ii));
    K3=dt*(ss_Ws.a*(ym(:,ii)+0.5*K2)+ss_Ws.b*r(ii));
    K4=dt*(ss_Ws.a*(ym(:,ii)+K3)+ss_Ws.b*r(ii));
    ym(:,ii+1)=ym(:,ii)+(K1+2*K2+2*K3+K4)/6;
    dym(:,ii)=(ym(:,ii+1)-ym(:,ii))/dt;
    
   % いろいろ必要な処理をしておく
   ym_true(ii)=ss_Ws.c*ym(:,ii); % yは状態量で，実際の出力をy_trueに設置
   dym_true(ii)=ss_Ws.c*dym(:,ii); % yは状態量で，実際の出力をy_trueに設置
   ddym_true(ii)=[1 0]*dym(:,ii);
   
   % omegaを求める
   omega1(:,ii)=[v1(:,ii);v2(:,ii);yf1(ii)];
   domega1(:,ii)=[dv1(:,ii);dv2(:,ii);dyf1(ii)];
   
   % z1を求める
   z1(ii)=y_true(ii)-ym_true(ii);
   
   % Ymを求めておく
   Ym(ii)=dym_true(ii)+lambda*ym_true(ii);
   dYm(ii)=ddym_true(ii)+lambda*dym_true(ii);
%    dz1(ii)=-lambda*z1(ii)+theta(:,ii)'*omega1(:,ii)+b0*uf(ii)-Ym(ii)+ep1(ii);
   
   % alpha1を求める
   fai(ii)=(c1-lambda)*z1(ii)+d1*z1(ii)+theta(:,ii)'*omega1(:,ii)-Ym(ii);
   alpha1(ii)=-p(ii)*fai(ii);
   
   % z2を求める
   z2(ii)=uf1(ii)-alpha1(ii);
   ganma1(ii)=-p(ii)*(c1-lambda+d1);
   ganmath1(:,ii)=-p(ii)*omega1(:,ii);
   aa=p(ii);
   bb=-((c1-lambda+d1)*z1(ii)+theta(:,ii)'*omega1(:,ii)-Ym(ii));
   cc=-p(ii)*theta(:,ii)';
   beta2(ii)=-lambda*uf1(ii)-ganma1(ii)*(-lambda*z1(ii)-Ym(ii))-aa*dYm(ii)-bb*g1*fai(ii)*z1(ii)-cc*domega1(:,ii);
%    dz2(ii)=u(ii)+beta2(ii)-ganma1(ii)*(theta(:,ii)'*omega1(:,ii)+b0**uf1(ii)+ep1(ii))-ganmath1(ii)*tauth2(ii)
   
   % alpha2を求める
   alpha2(ii)=-c2*z2(ii)-d2*ganma1(ii)^2*z2(ii)-beta2(ii)+ganma1(ii)*theta(:,ii)'*omega1(:,ii)+b0_hat(ii)*(ganma1(ii)*uf1(ii)-z1(ii))+ganmath1(:,ii)'*tauth2(:,ii);
   
   % 更新に使う変数
   tauth1(:,ii)=G*omega1(:,ii)*z1(ii);
   tauth2(:,ii)=tauth1(:,ii)-G*ganma1(ii)*omega1(:,ii)*z2(ii);
   taub2(ii)=g3*(z1(ii)-ganma1(ii)*uf1(ii))*z2(ii);
   
   % uを求める
%    u(ii)=-c2*z2(ii)-d2*ganma1(ii)^2*z2(ii)-beta2(ii)-b0_hat(ii)*z1(ii)+ganma1(ii)*theta(:,ii)'*omega1(:,ii)+b0_hat(ii)*ganma1(ii)*uf1(ii)+ganmath1(:,ii)'*tauth2(:,ii);
   u(ii)=alpha2(ii);

   % y(ii+1)を求める
    K1=dt*(ss_y.a*y(:,ii)+ss_y.b*u(ii));
    K2=dt*(ss_y.a*(y(:,ii)+0.5*K1)+ss_y.b*u(ii));
    K3=dt*(ss_y.a*(y(:,ii)+0.5*K2)+ss_y.b*u(ii));
    K4=dt*(ss_y.a*(y(:,ii)+K3)+ss_y.b*u(ii));
    y(:,ii+1)=y(:,ii)+(K1+2*K2+2*K3+K4)/6;
    
    % uf1を求める
    K1=dt*(-lambda*uf1(ii)+u(ii));
    K2=dt*(-lambda*(uf1(ii)+0.5*K1)+u(ii));
    K3=dt*(-lambda*(uf1(ii)+0.5*K2)+u(ii));
    K4=dt*(-lambda*(uf1(ii)+K3)+u(ii));
    uf1(ii+1)=uf1(ii)+(K1+2*K2+2*K3+K4)/6;
   
   % pを更新（オイラー）
   p(ii+1)=p(ii)+dt*g1*fai(ii)*z1(ii);
   
   % thetaを更新
   theta(:,ii+1)=theta(:,ii)+dt*tauth2(:,ii);
   
   % b0_hatを更新
   b0_hat(ii+1)=b0_hat(ii)+dt*taub2(ii);
   
end
toc
%% figure
figure('Name','相対次数が2の時のモデル規範型適応制御(バックステッピング)の応答例')
plot(t,ym_true(1:length(t)),t,y_true(1:length(t)),'--',t,z1,'-.','lineWidth',2);
legend('y_M','y','e')
ylim([-1 1])
grid on