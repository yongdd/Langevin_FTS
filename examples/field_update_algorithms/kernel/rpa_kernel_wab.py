import numpy as np
import matplotlib.pyplot as plt

# read data
data = np.load("fields_fts2d.npz")
print(data.files)

nx = data['nx'][0:2]
lx = data['lx'][0:2]
f = data['f']

w = data['w']
#w = [data['w_plus'] + data['w_minus'],data['w_plus'] - data['w_minus']]
w = np.reshape(w, (2, nx[0], nx[1]))

phi_a = np.reshape(data['phi_a'], nx)
phi_b = np.reshape(data['phi_b'], nx)

print(f)
print(phi_a.shape)

# arrays for debye functions
space_kx, space_ky = np.meshgrid(
    2*np.pi/lx[1]*np.arange(nx[1]//2+1),
    2*np.pi/lx[0]*np.concatenate([np.arange(nx[0]//2), nx[0]//2-np.arange(nx[0]//2)]))
mag2_k = (space_kx**2 + space_ky**2)
print(mag2_k.shape)

mag2_k[0,0] = 1.0e-5 # to prevent 'division by zero' error
g_aa_k = 2*(f*mag2_k+np.exp(-mag2_k*f)-1.0)/mag2_k**2
g_ab_k = (1.0-np.exp(-mag2_k*f))*(1.0-np.exp(-mag2_k*(1.0-f)))/mag2_k**2
g_bb_k = 2*((1.0-f)*mag2_k+np.exp(-mag2_k*(1.0-f))-1.0)/mag2_k**2
g_aa_k[0,0] = f**2
g_ab_k[0,0] = f*(1.0-f)
g_bb_k[0,0] = (1.0-f)**2

print(g_aa_k[0,0:5])
print(g_ab_k[0,0:5])
print(g_bb_k[0,0:5])

w_a_k = np.fft.rfftn(w[0])
w_b_k = np.fft.rfftn(w[1])
phi_1_k = -g_aa_k*w_a_k - g_ab_k*w_b_k
phi_2_k = -g_ab_k*w_a_k - g_bb_k*w_b_k

phi_1_k[0,0] = f*nx[0]*nx[1]
phi_2_k[0,0] = (1.0-f)*nx[0]*nx[1]

phi_1 = np.fft.irfftn(phi_1_k, nx)
phi_2 = np.fft.irfftn(phi_2_k, nx)


phi_rpa_k = np.zeros_like(phi_1_k)
#phi_rpa_k[0,0] = nx[0]*nx[1]
w_b_rpa_k = -(phi_rpa_k + w_a_k*(g_aa_k + g_ab_k))/(g_ab_k + g_bb_k)
w_b_rpa = np.fft.irfftn(w_b_rpa_k, nx)

print(np.mean(w_b_rpa_k), np.std(w_b_rpa_k))
print(np.mean(w_b_rpa), np.std(w_b_rpa))

print(np.mean(phi_a), np.std(phi_a))
print(np.mean(phi_1), np.std(phi_1))
print(np.mean(phi_b), np.std(phi_b))
print(np.mean(phi_2), np.std(phi_2))

print(np.std(phi_a-phi_1))
print(np.std(phi_b-phi_2))

fig, axes = plt.subplots(4,4, figsize=(20,20))
for i in range(0,4):
   for j in range(0,4): 
        axes[i,j].axis("off")

axes[0,0].imshow(w[0], vmin=-10.0, vmax=10.0, cmap="jet")
axes[0,1].imshow(w[1], vmin=-10.0, vmax=10.0, cmap="jet")
axes[0,2].imshow((w[0]+w[1])/2, vmin=-10.0, vmax=10.0, cmap="jet")
axes[0,3].imshow((w[0]-w[1])/2, vmin=-10.0, vmax=10.0, cmap="jet")
axes[1,1].imshow(w_b_rpa, vmin=-10.0, vmax=10.0, cmap="jet")

axes[2,0].imshow(phi_a, vmin=0.0, vmax=1.0, cmap="jet")
axes[2,1].imshow(phi_b, vmin=0.0, vmax=1.0, cmap="jet")
axes[2,2].imshow(phi_a+phi_b, vmin=0.5, vmax=1.5, cmap="jet")
axes[2,3].imshow(phi_a-phi_b, vmin=-1.0, vmax=1.0, cmap="jet")

axes[3,0].imshow(phi_1, vmin=0.0, vmax=1.0, cmap="jet")
axes[3,1].imshow(phi_2, vmin=0.0, vmax=1.0, cmap="jet")
axes[3,2].imshow(phi_1+phi_2, vmin=0.5, vmax=1.5, cmap="jet")
axes[3,3].imshow(phi_1-phi_2, vmin=-1.0, vmax=1.0, cmap="jet")


plt.subplots_adjust(left=0.01,bottom=0.01,
                    top=0.99,right=0.99,
                    wspace=0.01, hspace=0.01)
plt.savefig('kernel_field_wab.png')

