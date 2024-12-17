

title = "Galaxy Evolution Over Time"
filename = "files/tests/galaxy/milkyway_over_time.csv"


import imageio
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import os
import shutil

while not os.path.exists('files'):
    os.chdir('..')
print("Current working directory:", os.getcwd())



#################################
### PREPARE FILE ARCHITECTURE ###
#################################

folder = os.path.dirname(filename)
filename = os.path.basename(filename)
new_folder = filename.replace(".csv", "_gif")

if os.path.exists(os.path.join(folder, new_folder)):
    shutil.rmtree(os.path.join(folder, new_folder))

os.mkdir(os.path.join(folder, new_folder))

df = pd.read_csv(os.path.join(folder, filename))


#########################
### IMAGE PREPARATION ###
#########################



df = df[['x', 'y', 'z', 'phi']]
timeline = sorted(df['phi'].unique())

df_init = df[df['phi'] == 0]
x_min, x_max, y_min, y_max, z_min, z_max = df_init['x'].min(), df_init['x'].max(), df_init['y'].min(), df_init['y'].max(), df_init['z'].min(), df_init['z'].max()



from mpl_toolkits.mplot3d import Axes3D
from matplotlib.transforms import Affine2D

# Define the viewing parameters
R = max(x_max, y_max)  # Galaxy radius scale
h = 2*R              # Height for the viewpoint
viewpoints = 360       # Number of frames for a full rotation
angle_step = 360 / viewpoints

for i, t in enumerate(timeline):
    print(f"Rendering frame {i + 1}/{len(timeline)}", end="\r")
    df_t = df[df['phi'] == t]

    # Create a 3D figure
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Calculate current rotation angle
    theta = i * angle_step
    x_view = 2 * R * np.cos(np.radians(theta))
    y_view = 2 * R * np.sin(np.radians(theta))
    s = np.sin(np.radians(2 * theta)) # osciallate
    s /= (np.abs(s) + 0.001)**0.25 # make it go to 1 (slightly) faster
    z_view = h / 4 * (1+s)**0.5 # stay less time at bottom

    # Scatter the particles
    ax.scatter(df_t['x'], df_t['y'], df_t['z'], s=1, c='black', alpha=0.5)
    #ax.scatter(df_t['x'].iloc[0], df_t['y'].iloc[0], df_t['z'].iloc[0], s=20, c='red', zorder=10)  # Highlight SMBH

    # Set the camera/viewpoint
    ax.view_init(elev=np.degrees(np.arctan2(z_view, np.sqrt(x_view**2 + y_view**2))),
                 azim=np.degrees(np.arctan2(y_view, x_view)))

    # Axes limits
    ax.set_xlim(-R, R)
    ax.set_ylim(-R, R)
    ax.set_zlim(-R, R)
    ax.set_title(f"t={t}")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Save the frame
    plt.savefig(os.path.join(folder, new_folder, f"image_{i:03d}.png"))
    plt.close()

print("Progress... 100%")


##################
### CREATE GIF ###
##################


images = [f for f in os.listdir(os.path.join(folder, new_folder)) if f.endswith(".png")]
images = sorted(images, key=lambda x: (len(x), x))
gif_path = os.path.join(folder, filename.replace(".csv", ".gif"))

# Create a GIF from the PNG files
with imageio.get_writer(gif_path, mode='I', duration=0.03, loop=0) as writer:  # Adjust duration (in seconds) per frame as needed
    for i, png in enumerate(images):
        print(f"Progress... {i/len(images):.0%}", end="\r")
        image_path = os.path.join(folder, new_folder, png)
        writer.append_data(imageio.imread(image_path))

print(f"GIF created successfully at {gif_path}")

# let's clean up the folder
shutil.rmtree(os.path.join(folder, new_folder))