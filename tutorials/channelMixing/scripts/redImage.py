import cv2
import numpy as np

# Load the image
image = cv2.imread('../images/pitzDailiyPV_resizedEdited.tiff')

red = image.copy() 
red[np.sum(image, axis=2) > 500] = 0
red[:,:,0] = 0
red[:,:,1] = 0
gray = cv2.cvtColor(red, cv2.COLOR_BGR2GRAY)
ret,thresh = cv2.threshold(gray,10,255,0)

cv2.imwrite('../images/pitzDailiyPV_red.tiff', thresh)


blue = image.copy() 
blue[np.sum(image, axis=2) > 500] = 0
blue[:,:,1] = 0
blue[:,:,2] = 0
gray = cv2.cvtColor(blue, cv2.COLOR_BGR2GRAY)
ret,thresh = cv2.threshold(gray,10,255,0)

cv2.imwrite('../images/pitzDailiyPV_blue.tiff', thresh)
