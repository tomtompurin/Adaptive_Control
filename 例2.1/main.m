%% 適応制御　2章の例題（問題設定は2章と同じ）
clear
close all
%% 適応制御　1章の例題
clear
close all
%% 定数定義
% 規範モデル（設計変数）
am=1;
bm=1;
% 実際のモデル（未知です）
a=-0.5;
b=0.5;
% 適応側に関する定数（設計変数）
g1=5;
g2=5;
%% 入力生成
dt=0.01;
t=0:dt:50;
r=sin(t);
% r=sin(t)+3*sin(7*t); % 入力（指令値でも外乱でもよし）
%% 計算結果を入れる箱
ym=zeros(1,length(t)); % 規範モデルの応答
dym=zeros(1,length(t)); % 上の微分
y=zeros(1,length(t)); % 実際の制御対象の応答
dy=zeros(1,length(t)); % 上の微分
u=zeros(1,length(t)); % 制御対象への入力
theta1=zeros(1,length(t)); % フィードフォワードゲイン（1行目）とその微分値（2行目）
theta2=zeros(1,length(t)); % フィードバックゲイン（1行目）とその微分値（2行目）
e=zeros(1,length(t)); % 規範モデルの応答と実際の制御対象の応答の誤差
%% 計算（時間変化する各変数を計算）
% 適応則以外は変更がないことに注意
for ii=1:length(t)
    % ym(ii+1)をまず求める
    K1=dt*(bm*r(ii)-am*ym(ii));
    K2=dt*(bm*r(ii)-am*(ym(ii)+0.5*K1));
    K3=dt*(bm*r(ii)-am*(ym(ii)+0.5*K2));
    K4=dt*(bm*r(ii)-am*(ym(ii)+K3));
    ym(ii+1)=ym(ii)+(K1+2*K2+2*K3+K4)/6;
    dym(ii)=(ym(ii+1)-ym(ii))/dt;
    
    % 制御力uを求める
    % uやthetaの計算に制御対象の定数aやbが登場していないことがポイント
    u(ii)=theta1(ii)*r(ii)+theta2(ii)*y(ii);
    
    % エラーeを求める
    e(ii)=ym(ii)-y(ii);
    
    % y(ii+1)を求める
    K1=dt*(b*u(ii)-a*y(ii));
    K2=dt*(b*u(ii)-a*(y(ii)+0.5*K1));
    K3=dt*(b*u(ii)-a*(y(ii)+0.5*K2));
    K4=dt*(b*u(ii)-a*(y(ii)+K3));
    y(ii+1)=y(ii)+(K1+2*K2+2*K3+K4)/6;
    dy(ii)=(y(ii+1)-y(ii))/dt;
    
    % theta1(ii+1)を求める（P.27の適応則に基づく）
    theta1(ii+1)=theta1(ii)+dt*(g1*e(ii)*r(ii));
    
    % theta2(ii+1)を求める（P.27の適応則に基づく）
    theta2(ii+1)=theta2(ii)+dt*(g2*e(ii)*y(ii));
end
%% figure
figure('Name','安定論（リアプノフ法）に基づくモデル規範適応制御の応答例（定常ゲインと時定数が不明な1次系）')
plot(t,ym(1:length(t)),t,y(1:length(t)),'--',t,e,'-.','lineWidth',2);
legend('y_M','y','e')
grid on