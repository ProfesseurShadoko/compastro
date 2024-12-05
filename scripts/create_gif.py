

title = "Monopole Evolution"
filename = "files/output/evolve_monopole/mono_1000-iter.csv"


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

df = df[['x', 'y', 'phi']]
timeline = sorted(df['phi'].unique())

for i,t in enumerate(timeline):
    print(f"Progress... {i/len(timeline):.0%}", end="\r")
    df_t = df[df['phi'] == t]
    plt.figure(figsize=(10,10))
    plt.scatter(df_t['x'], df_t['y'])
    plt.xlim(-1, 1)
    plt.ylim(-1, 1)
    plt.title(f"t={t}")
    plt.suptitle(title)
    plt.savefig(os.path.join(folder, new_folder, f"image_{i}.png"))
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