```python
# https://github.com/worklifesg/Python-Computer-Vision-with-OpenCV-and-Deep-Learning
# downlaod the rainbow and crossword jpgs 
# make sure to go to terminal and install: pip install --upgrade pip & after that install: pip install opencv-python
```


```python
import cv2
import matplotlib.pyplot as plt 
%matplotlib inline
```


```python
img = cv2.imread('rainbow.jpg')
```


```python
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7fa463265550>




![png](output_3_1.png)



```python
# thresholding: cncel out certain background colors 
img = cv2.imread('rainbow.jpg', 0)
```


```python
plt.imshow(img, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7fa463204ed0>




![png](output_5_1.png)



```python
ret1, thresh1 = cv2.threshold(img, 127, 255, cv2.THRESH_BINARY)
```


```python
# there are 255 pixels in these images so 127 is nearly half the pixels without a decimal 
# everything below will be changed to one color and everything above will change to another color 
ret1
```




    127.0




```python
# everything that is lighter truns white, everything darker/gray turns black 
plt.imshow(thresh1, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7fa461175d50>




![png](output_8_1.png)



```python
# we can inverse the limits 
img2 = cv2.imread('rainbow.jpg', 0)
ret1, thresh1 = cv2.threshold(img2, 127, 255, cv2.THRESH_TRUNC)
plt.imshow(thresh1, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7fa4610e2590>




![png](output_9_1.png)



```python
# almost looks 3-D 
img3 = cv2.imread('rainbow.jpg', 0)
ret1, thresh1 = cv2.threshold(img3, 127, 255, cv2.THRESH_TOZERO)
plt.imshow(thresh1, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7fa4610c75d0>




![png](output_10_1.png)



```python
img_r = cv2.imread('crossword.jpg', 0)
plt.imshow(img_r, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7fa461025ed0>




![png](output_11_1.png)



```python
# creating this function to display our images 
def show_pic(img): 
    fig = plt.figure(figsize = (15,15))
    ax = fig.add_subplot(111)
    ax.imshow(img, cmap = 'gray')
```


```python
show_pic(img_r)
```


![png](output_13_0.png)



```python
# we cna cut out some extra steps by doing it this way: keep blac kand white, get rid of gray with this code 
ret, th1 = cv2.threshold(img_r, 127, 255, cv2.THRESH_BINARY)
show_pic(th1)
```


![png](output_14_0.png)



```python
ret, th1 = cv2.threshold(img_r, 200, 255, cv2.THRESH_BINARY)
show_pic(th1)
```


![png](output_15_0.png)



```python
# different technique to help get better results 
# in parenthases: (max value, adaptive method, threshold type, block size height, block size width)
th2 = cv2.adaptiveThreshold(img_r, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 11, 8)
show_pic(th2)
```


![png](output_16_0.png)



```python
# can blend to do different thresholds: essential layering 
blended = cv2.addWeighted(src1 = th1, alpha = 0.6, 
                          src2 = th2, beta = 0.4, gamma = 0)
show_pic(blended)
```


![png](output_17_0.png)



```python
# like the last one but got rid of some of the background 
th3 = cv2.adaptiveThreshold(img_r, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 11, 8)

blended = cv2.addWeighted(src1 = th1, alpha = 0.6, 
                          src2 = th3, beta = 0.4, gamma = 0)

show_pic(blended)
```


![png](output_18_0.png)

