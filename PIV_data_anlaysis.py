import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

Nx, Ny, Nz = 287, 190, 1000
Lx, Ly, d, u_sonic = 75.31583, 22.99752, 15, 315.35
all_u_data = np.zeros((Nx, Ny, Nz))
all_v_data = np.zeros((Nx, Ny, Nz))
u_bar = np.zeros((Nx, Ny))
v_bar = np.zeros((Nx, Ny))
x = np.linspace(0, Lx, Nx)
y = np.linspace(-Ly, Ly, Ny)
xnd = x/d
ynd = y/d
x_sonic_top = []
y_sonic_top = []
x_sonic_bottom = []
y_sonic_bottom = []
levels = np.linspace(0.0, 1.7, 18)
levels_rms = np.linspace(0.0, 0.5, 18)

def read_data(filename):
    data = pd.read_csv(filename)
    u_component = data['u [m/s]']
    v_component = data['v [m/s]']
    return u_component, v_component

def u_sample(N):
    u_bar = np.mean(all_u_data[:, :, :N], axis=2)
    v_bar = np.mean(all_v_data[:, :, :N], axis=2)
    return u_bar, v_bar

for k in range(1, 1000):
    filename = f'sample {k:04d} converted.csv'
    u, v = read_data(filename)
    all_u_data[:, :, k-1] = u.values.reshape(190, 287).T
    all_v_data[:, :, k-1] = v.values.reshape(190, 287).T

u_bar_100, v_bar_100 = u_sample(100)
u_bar_1000, v_bar_1000 = u_sample(1000)

u_rms_100 = np.sqrt(np.mean((all_u_data - u_bar_100[:, :, np.newaxis])**2, axis=2))
u_rms_1000 = np.sqrt(np.mean((all_u_data - u_bar_1000[:, :, np.newaxis])**2, axis=2))
v_rms_100 = np.sqrt(np.mean((all_v_data - v_bar_100[:, :, np.newaxis])**2, axis=2))
v_rms_1000 = np.sqrt(np.mean((all_v_data - v_bar_1000[:, :, np.newaxis])**2, axis=2))

for i in range(287):
    for j in range(95):
        if u_bar_1000[i, j]/u_sonic > 0.99 and u_bar_1000[i, j]/u_sonic < 1.01:
            x_sonic_top.append(xnd[i])
            y_sonic_top.append(ynd[j])

    for j in range(95, 190):
        if u_bar_1000[i, j]/u_sonic > 0.99 and u_bar_1000[i, j]/u_sonic < 1.01:
            x_sonic_bottom.append(xnd[i])
            y_sonic_bottom.append(ynd[j])

shock = np.argmax(u_rms_1000[200:, 95] / u_sonic) + 200

plt.figure(figsize=(10, 6))

plt.subplot(3, 2, 1)
plt.contourf(xnd, ynd, u_bar_100.T / u_sonic, levels=levels)
plt.colorbar()
plt.title('U Component Mean over 100 Samples')
plt.xlabel('x/d')
plt.ylabel('y/d')

plt.subplot(3, 2, 2)
plt.contourf(xnd, ynd, u_bar_1000.T / u_sonic, levels=levels)
plt.colorbar()
plt.title('U Component Mean over 1000 Samples')
plt.xlabel('x/d')
plt.ylabel('y/d')

plt.subplot(3, 2, 3)
plt.contourf(xnd, ynd, u_rms_100.T / u_sonic, levels=levels_rms)
plt.colorbar()
plt.title('U rms over 100 Samples')
plt.xlabel('x/d')
plt.ylabel('y/d')

plt.subplot(3, 2, 4)
plt.contourf(xnd, ynd, u_rms_1000.T / u_sonic, levels=levels_rms)
plt.colorbar()
plt.title('U rms over 1000 Samples')
plt.xlabel('x/d')
plt.ylabel('y/d')

plt.subplot(3, 2, 5)
plt.contourf(xnd, ynd, v_rms_100.T / u_sonic, levels=levels_rms)
plt.colorbar()
plt.title('V rms over 100 Samples')
plt.xlabel('x/d')
plt.ylabel('y/d')

plt.subplot(3, 2, 6)
plt.contourf(xnd, ynd, v_rms_1000.T / u_sonic, levels=levels_rms)
plt.colorbar()
plt.title('V rms over 1000 Samples')
plt.xlabel('x/d')
plt.ylabel('y/d')

plt.savefig('velocity statistics.png')
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.plot(xnd, u_rms_100[:, 95] / u_sonic, label='U rms at y/d=0')
plt.plot(xnd, u_rms_100[:, 62] / u_sonic, label='U rms at y/d=0.5')
plt.title('U rms vs x (100 samples)')
plt.xlabel('x/d')
plt.ylabel('U rms / U sonic')
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(xnd, u_rms_1000[:, 95] / u_sonic, label='U rms at y/d=0')
plt.plot(xnd, u_rms_1000[:, 62] / u_sonic, label='U rms at y/d=0.5')
plt.title('U rms vs x (1000 samples)')
plt.xlabel('x/d')
plt.ylabel('U rms / U sonic')
plt.legend()

plt.savefig('u mean profiles.png')
plt.tight_layout()
plt.show()

plt.figure(figsize=(14, 5))

plt.subplot(1, 2, 1)
plt.contourf(xnd, ynd, u_rms_1000.T / u_sonic, levels=levels)
plt.plot(x_sonic_top, y_sonic_top, color='red', label='Sonic Curve')
plt.plot(x_sonic_bottom, y_sonic_bottom, color='red')
plt.axvline(xnd[shock], linestyle='--', label='Shock Location')
plt.colorbar()
plt.title('Sonic Curve in Jet Shear Layer')
plt.xlabel('x/d')
plt.ylabel('y/d')
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(xnd, u_rms_1000[:, 95] / u_sonic, label='U rms at y/d=0')
plt.axvline(xnd[shock], linestyle='--', label='Shock Location')
plt.title('U rms vs x (1000 samples)')
plt.xlabel('x/d')
plt.ylabel('U rms / U sonic')
plt.legend()

plt.savefig('sonic curve and shock location.png')
plt.tight_layout()
plt.show()
