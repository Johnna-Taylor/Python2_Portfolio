# Python2_Portfolio
This is the portfolio of python that I learned during BISC450C-V84

## Open CV 
### Part 3

In this class section we learned how to do gray scale (black/white) with thresholds 

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




![output_3_1](https://github.com/user-attachments/assets/ea017b69-5696-46f3-a3c0-82ea61632211)




```python
# thresholding: cncel out certain background colors 
img = cv2.imread('rainbow.jpg', 0)
```


```python
plt.imshow(img, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7fa463204ed0>



![output_5_1](https://github.com/user-attachments/assets/aef17f19-9c7e-470a-9df3-97314259bc9a)





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




![output_8_1](https://github.com/user-attachments/assets/3b035f2c-8cd8-42a7-b3c7-9f424268c801)




```python
# we can inverse the limits 
img2 = cv2.imread('rainbow.jpg', 0)
ret1, thresh1 = cv2.threshold(img2, 127, 255, cv2.THRESH_TRUNC)
plt.imshow(thresh1, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7fa4610e2590>




![output_9_1](https://github.com/user-attachments/assets/2c5affd3-4fc0-4dbd-9204-5ccc90ebc09d)




```python
# almost looks 3-D 
img3 = cv2.imread('rainbow.jpg', 0)
ret1, thresh1 = cv2.threshold(img3, 127, 255, cv2.THRESH_TOZERO)
plt.imshow(thresh1, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7fa4610c75d0>




![output_10_1](https://github.com/user-attachments/assets/8be4f704-e5cb-4bf8-b890-021ec865c437)




```python
img_r = cv2.imread('crossword.jpg', 0)
plt.imshow(img_r, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7fa461025ed0>




![output_11_1](https://github.com/user-attachments/assets/0def881d-2516-4c1e-97be-54ad6a1686f9)




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


![output_13_0](https://github.com/user-attachments/assets/adce4b78-046a-4a03-9f77-0431f0fb5792)




```python
# we cna cut out some extra steps by doing it this way: keep blac kand white, get rid of gray with this code 
ret, th1 = cv2.threshold(img_r, 127, 255, cv2.THRESH_BINARY)
show_pic(th1)
```


![output_14_0](https://github.com/user-attachments/assets/23b8cdea-67fa-4f22-b2ca-f5abb3efa550)




```python
ret, th1 = cv2.threshold(img_r, 200, 255, cv2.THRESH_BINARY)
show_pic(th1)
```


![output_15_0](https://github.com/user-attachments/assets/0a138112-dcd7-462a-8657-d419405c682e)




```python
# different technique to help get better results 
# in parenthases: (max value, adaptive method, threshold type, block size height, block size width)
th2 = cv2.adaptiveThreshold(img_r, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 11, 8)
show_pic(th2)
```


![output_16_0](https://github.com/user-attachments/assets/87e3b2a2-6515-4714-ab84-1f9d11db5a03)




```python
# can blend to do different thresholds: essential layering 
blended = cv2.addWeighted(src1 = th1, alpha = 0.6, 
                          src2 = th2, beta = 0.4, gamma = 0)
show_pic(blended)
```


![output_17_0](https://github.com/user-attachments/assets/0fbfb928-8d13-4467-aead-6640aca0de1d)




```python
# like the last one but got rid of some of the background 
th3 = cv2.adaptiveThreshold(img_r, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 11, 8)

blended = cv2.addWeighted(src1 = th1, alpha = 0.6, 
                          src2 = th3, beta = 0.4, gamma = 0)

show_pic(blended)
```


![output_18_0](https://github.com/user-attachments/assets/d1849740-a31e-45e3-afaa-507d4178e0ce)


