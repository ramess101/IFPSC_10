import matplotlib.pyplot as plt, pandas, numpy as np

df = pandas.read_excel('IOCTANE.xlsx')
df_fit = df[(df['-sr/R'] > 3.5) & (df.used) & (df['T'] > 240.0)].copy()
df = df[(df['-sr/R'] > 3.5) & (df.used)].copy()
df_room = df[(df['T'] > 283.0) & (df['T'] < 299.0)].copy()
df_low = df[(df['T'] < 240.0)].copy()
df_high = df[(df['T'] > 299.0)].copy()

cc = np.polyfit(df_fit['-sr/R'], np.log10(df_fit['eta/eta0']), 3)
etarfit = 10**(np.polyval(cc, df['-sr/R']))
etarfit_room = 10**(np.polyval(cc, df_room['-sr/R']))
etarfit_low = 10**(np.polyval(cc, df_low['-sr/R']))
etarfit_high = 10**(np.polyval(cc, df_high['-sr/R']))
xx = np.linspace(np.min(df['-sr/R']), np.max(df['-sr/R'])+1,1000)

dev = np.abs(etarfit/df['eta/eta0']-1)*100
dev_room = np.abs(etarfit_room/df_room['eta/eta0']-1)*100
dev_low = np.abs(etarfit_low/df_low['eta/eta0']-1)*100
dev_high = np.abs(etarfit_high/df_high['eta/eta0']-1)*100

fig,(ax1,ax2,ax3) = plt.subplots(3,1,figsize=(6,12),sharex=True)

#ax1.plot(df['-sr/R'], df['eta/eta0'],'o')
ax1.plot(df_low['-sr/R'], df_low['eta/eta0'],'bs',alpha=0.5)
ax1.plot(df_high['-sr/R'], df_high['eta/eta0'],'r^',alpha=0.5)
ax1.plot(df_room['-sr/R'], df_room['eta/eta0'],'go',alpha=0.5)
ax1.set_ylabel(r'$\eta/\eta_{\rho\to 0}$')

#ax2.plot(df['-sr/R'], df['eta/eta0'],'o')
ax2.plot(df_low['-sr/R'], df_low['eta/eta0'],'bs',alpha=0.5)
ax2.plot(df_high['-sr/R'], df_high['eta/eta0'],'r^',alpha=0.5)
ax2.plot(df_room['-sr/R'], df_room['eta/eta0'],'go',alpha=0.5)
ax2.plot(xx, 10**(np.polyval(cc, xx)),'k--')
ax2.set_yscale('log')
ax2.set_ylabel(r'$\eta/\eta_{\rho\to 0}$')

#ax3.plot(df['-sr/R'],dev,'o')
ax3.plot(df_low['-sr/R'],dev_low,'bs',alpha=0.5)
ax3.plot(df_high['-sr/R'],dev_high,'r^',alpha=0.5)
ax3.plot(df_room['-sr/R'],dev_room,'go',alpha=0.5)
ax3.set_xlabel(r'$-s^{\rm r}/R$')
ax3.set_ylabel(r'$|(\eta/\eta_{\rm fit}-1)\times 100|$')
plt.legend(['T < 240 K','T > 298 K','283 < T/K < 298'])
plt.tight_layout(pad=0.2)
plt.savefig('IOCTANE.pdf')
plt.show()

fig,(ax1,ax2,ax3) = plt.subplots(3,1,figsize=(6,12),sharex=True)

p = ax1.scatter(df['-sr/R'], df['eta/eta0'],c=df['T'],cmap='rainbow')
ax1.set_ylabel(r'$\eta/\eta_{\rho\to 0}$')

p = ax2.scatter(df['-sr/R'], df['eta/eta0'],c=df['T'],cmap='rainbow')
ax2.plot(xx, 10**(np.polyval(cc, xx)))
ax2.set_yscale('log')
ax2.set_ylabel(r'$\eta/\eta_{\rho\to 0}$')

dev = np.abs(etarfit/df['eta/eta0']-1)*100
dev_room = np.abs(etarfit_room/df_room['eta/eta0']-1)*100

p = ax3.scatter(df['-sr/R'],dev,c=df['T'],cmap='rainbow')
ax3.set_xlabel(r'$-s^{\rm r}/R$')
ax3.set_ylabel(r'$|(\eta/\eta_{\rm fit}-1)\times 100|$')
plt.tight_layout(pad=0.6)

cb = plt.colorbar(p,ax=ax1,pad=0.01)
cb.set_label('Temperature (K)')

plt.savefig('IOCTANE_2.pdf')
plt.show()