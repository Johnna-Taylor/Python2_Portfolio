```python
# first go to terminal and install: pip install --upgrade pip
# make sure to install: pip install opencv-python
import cv2
import numpy as np
```


```python
import matplotlib.pyplot as plt 
%matplotlib inline
```


```python
img = cv2.imread("mushroom.jpg")
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7f2652952650>




![png](output_2_1.png)



```python
# setting the color ratio for the output pictures; 25 is the max but we want two colors 
edges = cv2.Canny(image =img, threshold1 = 127, threshold2 = 127)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7f2652083d90>




![png](output_3_1.png)



```python
# this gives us our median color value, this tells us that our image is darker 
med_value = np.median(img)
med_value
```




    58.0




```python
# creating new values color variables
lower = int(max(0,0.7*med_value))
upper = int(min(25,1.3*med_value))

edges = cv2.Canny(img, threshold1 = lower, threshold2 = upper)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7f2652001fd0>




![png](output_5_1.png)



```python
edges = cv2.Canny(image = img, threshold1 = lower, threshold2 = upper + 100)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7f2651f6fd50>




![png](output_6_1.png)



```python
# blurring as a solution; creates lower resolution with pixels  
blurred_img = cv2.blur(img, ksize = (5,5))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower, 
                 threshold2 = upper)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7f2651ee2550>




![png](output_7_1.png)



```python
blurred_img = cv2.blur(img, ksize = (7,7))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower, 
                 threshold2 = upper)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7f2651ec98d0>




![png](output_8_1.png)



```python
# can increase the upper threshold limit 
blurred_img = cv2.blur(img, ksize = (5,5))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower, 
                 threshold2 = upper + 50)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7f2651e3a990>




![png](output_9_1.png)



```python
blurred_img = cv2.blur(img, ksize = (5,5))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower, 
                 threshold2 = upper + 100)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7f2651da7a90>




![png](output_10_1.png)



```python
blurred_img = cv2.blur(img, ksize = (7,7))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower, 
                 threshold2 = upper + 60)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7f2651d16bd0>




![png](output_11_1.png)



```python
blurred_img = cv2.blur(img, ksize = (7,7))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower, 
                 threshold2 = upper + 60)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7f2651d03b50>




![png](output_12_1.png)



```python
# I was just messing around trying to see how well i could get mine to outline
blurred_img = cv2.blur(img, ksize = (1,1))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower, 
                 threshold2 = upper + 145)
                  
plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7f2651c6bc90>




![png](output_13_1.png)

