import pyMRAW
import cv2
import numpy as np
import matplotlib.pyplot as plt


# This loads all frames into RAM
images, info = pyMRAW.load_video(r"A:\OneDrive - Monash University\Uni\HPR\FYP\Testing data\Sample\test_footage.cih")
print("Loaded video shape:", images.shape)
print("Info:", info)
print("Parsed header:")
for k, v in info.items():
    print(f"{k}: {v}")

first_image = images[0]

plt.figure()
plt.imshow(first_image, cmap='gray')
plt.show()
