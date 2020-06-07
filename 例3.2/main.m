%% 適応制御　3章の例題
% 相対次数が2の時（3次以降もこのプログラムを少し変更すればできます）
clear
close all
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

%% zeta1(適応則更新用の変数)を求めるのに使う状態空間表現
Az=[zeros(5,5) eye(5);
    -eye(5) -2*eye(5)];
Bz=[zeros(5,5);eye(5)];
Cz=[eye(5) zeros(5,5)];
%% その他必要な定数
b0=ss_y.c*ss_y.A^(nn-1)*ss_y.b;
lambda=1; % >0
g1=1; % >0

%% 状態変数フィルタ（設計変数）
% (F,g)が可制御になるように選ぶ
F=[-1 0;
    0 -2];
g=[1;1];

%% 適応ゲイン 今回はn次元のフィルタを用いた
G=eye(5);

%% 入力
dt=0.2; % 刻み時間はマシンスペックとか計算誤差を見ながら決めてます。
t=0:dt:40;
r=sin(t);

%% 計算結果を入れる箱
ym=zeros(2,length(t)); % 規範モデルの応答
ym_true=zeros(1,length(t)); % 規範モデルの応答
y=zeros(2,length(t)); % 実際の制御対象の応答
y_true=zeros(1,length(t)); % 実際の制御対象の応答
u=zeros(1,length(t)); % 制御対象への入力
omega=zeros(2*n+1,length(t)); % 適応や入力に使う変数群
e=zeros(1,length(t)); % 規範モデルの応答と実際の制御対象の応答の誤差
e_hat=zeros(1,length(t)); % 規範モデルの応答と実際の制御対象の推定誤差
ea=zeros(1,length(t)); % 規範モデルの応答と実際の制御対象の拡張誤差
theta=zeros(2*n+1,length(t)); % フィードバックゲイン
v1=zeros(2,length(t)); % 状態変数フィルタの状態量
v2=zeros(2,length(t)); % 状態変数フィルタの状態量
b0_hat=zeros(1,length(t)); % b0の推定値（b0を使えればいいが，制御対象は未知なので）
zeta1=zeros(2*(2*n+1),length(t)); % 適応則に使う中間変数
zeta_true1=zeros((2*n+1),length(t)); % 適応則に使う中間変数
zeta2=zeros(2,length(t)); % 適応則に使う中間変数
zeta_true2=zeros(1,length(t)); % 適応則に使う中間変数
m=zeros(1,length(t)); % 適応則に使う中間変数
z=zeros(1,length(t)); % 適応則に使う中間変数
%% 計算 時間変化していく各変数を計算
tic
for ii=1:length(t)
    % e(ii)を計算しておく
    y_true(ii)=ss_y.c*y(:,ii); % yは状態量で，実際の出力をy_trueに設置
    ym_true(ii)=ss_Ws.c*ym(:,ii); % yは状態量で，実際の出力をy_trueに設置
    e(ii)=ym_true(ii)-y_true(ii);
    
    % e_hat(ii)を計算しておく
    e_hat(ii)=b0_hat(ii)*z(ii);
    
    % ea(ii)を計算しておく
    ea(ii)=e(ii)-e_hat(ii);
    
    % omega(ii)を計算しておく
    omega(:,ii)=[r(ii);v1(:,ii);v2(:,ii)];
    
    % uを計算しておく
    u(ii)=theta(:,ii)'*omega(:,ii);
    
    % y(ii+1)を求める
    K1=dt*(ss_y.a*y(:,ii)+ss_y.b*u(ii));
    K2=dt*(ss_y.a*(y(:,ii)+0.5*K1)+ss_y.b*u(ii));
    K3=dt*(ss_y.a*(y(:,ii)+0.5*K2)+ss_y.b*u(ii));
    K4=dt*(ss_y.a*(y(:,ii)+K3)+ss_y.b*u(ii));
    y(:,ii+1)=y(:,ii)+(K1+2*K2+2*K3+K4)/6;
    
    % v1(ii+1)を求める
    K1=dt*(g*u(ii)+F*v1(:,ii));
    K2=dt*(g*u(ii)+F*(v1(:,ii)+0.5*K1));
    K3=dt*(g*u(ii)+F*(v1(:,ii)+0.5*K2));
    K4=dt*(g*u(ii)+F*(v1(:,ii)+K3));
    v1(:,ii+1)=v1(:,ii)+(K1+2*K2+2*K3+K4)/6;
    
    % v2(ii+1)を求める
    K1=dt*(g*y_true(ii)+F*v2(:,ii));
    K2=dt*(g*y_true(ii)+F*(v2(:,ii)+0.5*K1));
    K3=dt*(g*y_true(ii)+F*(v2(:,ii)+0.5*K2));
    K4=dt*(g*y_true(ii)+F*(v2(:,ii)+K3));
    v2(:,ii+1)=v2(:,ii)+(K1+2*K2+2*K3+K4)/6;
    
    % zeta_true1(ii+1)を求める
    K1=dt*(Az*zeta1(:,ii)+Bz*omega(:,ii));
    K2=dt*(Az*(zeta1(:,ii)+0.5*K1)+Bz*omega(:,ii));
    K3=dt*(Az*(zeta1(:,ii)+0.5*K2)+Bz*omega(:,ii));
    K4=dt*(Az*(zeta1(:,ii)+K3)+Bz*omega(:,ii));
    zeta1(:,ii+1)=zeta1(:,ii)+(K1+2*K2+2*K3+K4)/6;
    zeta_true1(:,ii+1)=Cz*zeta1(:,ii+1);
    
    % zeta_true1(ii+1)を求める
    K1=dt*(ss_Ws.a*zeta2(:,ii)+ss_Ws.b*u(ii));
    K2=dt*(ss_Ws.a*(zeta2(:,ii)+0.5*K1)+ss_Ws.b*u(ii));
    K3=dt*(ss_Ws.a*(zeta2(:,ii)+0.5*K2)+ss_Ws.b*u(ii));
    K4=dt*(ss_Ws.a*(zeta2(:,ii)+K3)+ss_Ws.b*u(ii));
    zeta2(:,ii+1)=zeta2(:,ii)+(K1+2*K2+2*K3+K4)/6;
    zeta_true2(:,ii+1)=ss_Ws.c*zeta2(:,ii+1);
    
    % m(ii)を求める
    m(ii)=lambda+zeta_true1(:,ii+1)'*zeta_true1(:,ii+1);
    
    % z(ii)を求める
    z(ii)=theta(:,ii)'*zeta_true1(:,ii)-zeta_true2(:,ii);
    
    % b0_hat(ii+1)を求める（オイラー法）
    b0_hat(ii+1)=b0_hat(ii)+dt*g1*z(ii)/m(ii)*ea(ii);
    
    % theta(ii+1)を求める（オイラー法）
    theta(:,ii+1)=theta(:,ii)+dt*sign(b0)*G*zeta_true1(:,ii)/m(ii)*ea(ii);
    
    % ym(ii+1)を求める
    K1=dt*(ss_Ws.a*ym(:,ii)+ss_Ws.b*r(ii));
    K2=dt*(ss_Ws.a*(ym(:,ii)+0.5*K1)+ss_Ws.b*r(ii));
    K3=dt*(ss_Ws.a*(ym(:,ii)+0.5*K2)+ss_Ws.b*r(ii));
    K4=dt*(ss_Ws.a*(ym(:,ii)+K3)+ss_Ws.b*r(ii));
    ym(:,ii+1)=ym(:,ii)+(K1+2*K2+2*K3+K4)/6;
end
toc
%% figure
figure('Name','相対次数が2の時のモデル規範型適応制御の応答例')
plot(t,ym_true(1:length(t)),t,y_true(1:length(t)),'--',t,e,'-.','lineWidth',2);
legend('y_M','y','e')
grid on