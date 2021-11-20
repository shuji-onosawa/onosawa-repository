import numpy as np
import matplotlib.pyplot as plt
plt.figure()


LX, LY=10,10   # メッシュのためのパラメータ
gridwidth=0.1
X, Y= np.meshgrid(np.arange(-LX, LX, gridwidth), np.arange(-LY, LY,gridwidth))  # メッシュ生成

#物理パラメータ
beta = 1*10**2
alpha = 9*10**4
Re = 6.4*10**5
Beq = 3*10**(-5)
R = np.sqrt(X**2+Y**2) # 原点からの距離

# ベクトル関数 F(U(x,y), V(x,y))を定義。静電場の表式を用いる。
B = Beq*(Re/R)**3

Ex = beta*(R + Y**2/R)/(B*Re**2) - Re*alpha*Y/(B*R**3)
Ey = -beta*X*Y/(B*R*Re**2) + Re*alpha*X/(B*R**3)

plt.streamplot(X,Y,Ex,Ey) # 流線をプロット

plt.xlim([-LX,LX]) # 描くXの範囲
plt.ylim([-LY,LY]) # 描くyの範囲

# グラフ描画
plt.grid()
plt.draw()
plt.show()

