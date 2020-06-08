%% 適応制御　4章の例題
% 例4.2 切り替え型σ修正法
clear
close all
%% システム定義
%% 定数定義
% 規範モデル（設計変数）
am=1;
bm=1;
% 実際のモデル（未知です）
a=-1;
b=1;
% 適応側に関する定数（設計変数）
g1=1;

% σ修正法のσ
sigma0=1;

% 切り替え則に使う定数
M0=5;
%% 入力生成
dt=0.01;
t=0:dt:50;
r=zeros(1,length(t));
d=(1+t).^(-1/5).*(5-(1+t).^(-1/5)-0.4*(1+t).^(-6/5)); % 図4.6の場合の外乱
% d=zeros(1,length(t)); % 図4.7の場合の外乱
%% 計算結果を入れる箱
ym=zeros(1,length(t)); % 規範モデルの応答
y=zeros(1,length(t)); % 実際の制御対象の応答
u=zeros(1,length(t)); % 制御対象への入力
theta=zeros(1,length(t)); % 適応パラ
e=zeros(1,length(t)); % 規範モデルの応答と実際の制御対象の応答の誤差
sigma=zeros(1,length(t)); % 切り替え型のσ

%% 各パラメータの初期値
y(1)=1;
theta(1)=-5;

%% 計算（時間変化する各変数を計算）
for ii=1:length(t)
    % ym(ii+1)をまず求める
    K1=dt*(bm*r(ii)-am*ym(ii));
    K2=dt*(bm*r(ii)-am*(ym(ii)+0.5*K1));
    K3=dt*(bm*r(ii)-am*(ym(ii)+0.5*K2));
    K4=dt*(bm*r(ii)-am*(ym(ii)+K3));
    ym(ii+1)=ym(ii)+(K1+2*K2+2*K3+K4)/6;
    
    % 制御力uを求める
    % uやthetaの計算に制御対象の定数aやbが登場していないことがポイント
    u(ii)=r(ii)+theta(ii)*y(ii);
    
    % エラーeを求める
    e(ii)=ym(ii)-y(ii);
    
    % y(ii+1)を求める
    K1=dt*(-a*y(ii)+b*u(ii)+d(ii));
    K2=dt*(-a*(y(ii)+0.5*K1)+b*u(ii)+d(ii));
    K3=dt*(-a*(y(ii)+0.5*K2)+b*u(ii)+d(ii));
    K4=dt*(-a*(y(ii)+K3)+b*u(ii)+d(ii));
    y(ii+1)=y(ii)+(K1+2*K2+2*K3+K4)/6;
    
    % 2段階切り替え（コメントアウトで2段階切り替えと3段階切り替えを変更可能）
    if abs(theta(ii))>M0
        sigma(ii)=sigma0;
    else
        sigma(ii)=0;
    end
    
%     3段階切り替え（コメントアウトで2段階切り替えと3段階切り替えを変更可能）
%     if abs(theta(ii))>2*M0
%         sigma(ii)=sigma0;
%     elseif abs(theta(ii))<M0
%         sigma(ii)=0;
%     else
%         sigma(ii)=sigma0*(abs(theta(ii))/M0-1);
%     end
    
    % theta(ii+1)を求める（例4.1との差分はこれ）
    K1=dt*(-sigma(ii)*theta(ii)+g1*y(ii)*e(ii));
    K2=dt*(-sigma(ii)*(theta(ii)+0.5*K1)+g1*y(ii)*e(ii));
    K3=dt*(-sigma(ii)*(theta(ii)+0.5*K2)+g1*y(ii)*e(ii));
    K4=dt*(-sigma(ii)*(theta(ii)+K3)+g1*y(ii)*e(ii));
    theta(ii+1)=theta(ii)+(K1+2*K2+2*K3+K4)/6;
end
%% figure
figure('Name','切り替え型σ修正法の場合の応答')
plot(t,y(1:length(t)),t,theta(1:length(t)),'-.','lineWidth',2);
legend('y','θ')
xlabel('Time (s)')
ylabel('応答')
grid on