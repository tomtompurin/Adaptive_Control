%% 適応制御　3章の例題
% 相対次数が1次の場合（特殊な例）
clear
close all
%% システム定義
s=tf('s');
tf_y=(s+1)/s^2; % 制御対象のモデル
tf_ym=1/(s+1); % 規範モデル
Ws=s+1; % 規範モデルの分母
nn=1; % システムの相対次数
n=2; % システムの次数
%% 制御対象の実現
ss_y=ss(tf_y);

%% その他必要な定数
b0=ss_y.c*ss_y.A^(nn-1)*ss_y.b;

%% 状態変数フィルタ（設計変数）
F=[-1 0;
    0 -2];
g=[1;1];

%% 適応ゲイン 今回はn次元のフィルタを用いた
G=100*eye(5);

%% 入力
dt=0.01;
t=0:dt:20;
r=sin(t);

%% 計算結果を入れる箱
ym=zeros(1,length(t)); % 規範モデルの応答
y=zeros(2,length(t)); % 実際の制御対象の応答
y_true=zeros(1,length(t)); % 実際の制御対象の応答
u=zeros(1,length(t)); % 制御対象への入力
omega=zeros(5,length(t)); % 適応や入力に使う変数群
e=zeros(1,length(t)); % 規範モデルの応答と実際の制御対象の応答の誤差
theta=zeros(5,length(t)); % フィードバックゲイン
v1=zeros(2,length(t)); % 状態変数フィルタの状態量
v2=zeros(2,length(t)); % 状態変数フィルタの状態量

%% 計算 時間変化していく各定数を計算
for ii=1:length(t)  
    % e(ii)を計算しておく
    y_true(ii)=ss_y.c*y(:,ii); % yは状態量で，実際の出力をy_trueに設置
    e(ii)=ym(ii)-y_true(ii);
    
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
    
    % theta(ii+1)を求める(オイラー法)
    theta(:,ii+1)=theta(:,ii)+dt*sign(b0)*G*omega(:,ii)*e(ii);
    
    % ym(ii+1)を求める
    K1=dt*(r(ii)-ym(ii));
    K2=dt*(r(ii)-(ym(ii)+0.5*K1));
    K3=dt*(r(ii)-(ym(ii)+0.5*K2));
    K4=dt*(r(ii)-(ym(ii)+K3));
    ym(ii+1)=ym(ii)+(K1+2*K2+2*K3+K4)/6;
end

%% figure
figure('Name','相対次数が1の時のモデル規範型適応制御の応答例')
plot(t,ym(1:length(t)),t,y_true(1:length(t)),'--',t,e,'-.','lineWidth',2);
legend('y_M','y','e')
grid on