import cv2

# Load the image
image = cv2.imread('../images/pitzDailiyPV_croped.png')

# Define the desired dimensions for the resized image
width = 984
height = 160

# Resize the image
resized_image = cv2.resize(image, (width, height))
gray = cv2.cvtColor(resized_image, cv2.COLOR_BGR2GRAY)
ret,thresh = cv2.threshold(gray,100,255,0)

# Save the resized image
cv2.imwrite('../images/pitzDailiyPV_resized.tiff', thresh)
