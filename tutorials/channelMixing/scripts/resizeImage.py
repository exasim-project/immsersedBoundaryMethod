import cv2
import numpy as np
import os

if not os.path.isfile('../images/fge.tiff'):
    cv2.imwrite('../images/fge.tiff', np.zeros([160, 400, 3], dtype=np.uint8))
    exit()

# Load the image
image = cv2.imread('../images/fge.tiff')

# Define the desired dimensions for the resized image
width = 400
height = 160

# Resize the image
resized_image = cv2.resize(image, (width, height))
gray = cv2.cvtColor(resized_image, cv2.COLOR_BGR2GRAY)
ret,thresh = cv2.threshold(gray,100,255,0)

# Save the resized image
cv2.imwrite('../images/fge_resized.tiff', thresh)
