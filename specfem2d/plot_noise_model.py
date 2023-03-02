import matplotlib.pyplot as plt
import cmasher as cmr
import numpy as np
import os


def _read(filename, dtype='float32'):
    nbytes = os.path.getsize(filename)
    with open(filename, 'rb') as file:
        # read size of record
        file.seek(0)
        n = np.fromfile(file, dtype='int32', count=1)[0]

        if n == nbytes-8:
            file.seek(4)
            data = np.fromfile(file, dtype=dtype)
            return data[:-1]
        else:
            file.seek(0)
            data = np.fromfile(file, dtype=dtype)
            return data


# parameters
coor_dir = './sim_reference/OUTPUT_FILES'
model_dir = './sim_reference/NOISE_TOMOGRAPHY'
nproc = 16

plt_fsize = 16
plt_title = None
figname = None

# DO NOT EDIT BELOW THIS LINE
# ===========================
xcoor = np.array([])
zcoor = np.array([])
model = np.array([])

for i in range(0, nproc):
    mesh_glob = np.loadtxt(os.path.join(coor_dir, f'mesh_glob{i:06}'))

    xcoor_proc = mesh_glob[:, 0]
    zcoor_proc = mesh_glob[:, 1]

    model_proc = _read(os.path.join(model_dir, f'proc{i:06}_mask_noise.bin'))

    xcoor = np.append(xcoor, xcoor_proc)
    zcoor = np.append(zcoor, zcoor_proc)
    model = np.append(model, model_proc)

# convert to km
xcoor /= 1000.0
zcoor /= 1000.0

model /= np.max(model)

# figure
fig, ax = plt.subplots()

cmap = cmr.cosmic
im = ax.tripcolor(xcoor, zcoor, model, cmap=cmap, linewidth=0.0,
                  edgecolor='none')

cbar = plt.colorbar(mappable=im, pad=0.01, fraction=0.03)
cbar.ax.tick_params(labelsize=plt_fsize)
cbar.set_label('Normalized PSD', labelpad=10, fontsize=plt_fsize)

if plt_title:
    plt.title(plt_title, fontsize=plt_fsize)

plt.xlabel('X [km]', fontsize=plt_fsize)
plt.ylabel('Y [km]', fontsize=plt_fsize)
plt.xticks(fontsize=plt_fsize*0.8)
plt.yticks(fontsize=plt_fsize*0.8)
plt.xlim(xmin=xcoor.min(), xmax=xcoor.max())
plt.ylim(ymin=zcoor.min(), ymax=zcoor.max())

if not figname:
    plt.show()
else:
    plt.savefig(figname, dpi=600, bbox_inches='tight', facecolor='white')

plt.show()
