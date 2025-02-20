```python
# make sure to go to terminal and install: pip install --upgrade pip 
# & then install: pip install opencv-python
```


```python
import cv2 
import numpy as np
import matplotlib.pyplot as plt 
%matplotlib inline 
```


```python
# upload a picture of a green & white chessboard with no border or peices 
# upload a picture of a chessboard with a border and pieces 
```


```python
# how to upload the green chessboard we have saved 
flat_chess = cv2.imread('chessboard-green.png')
flat_chess = cv2.cvtColor(flat_chess, cv2.COLOR_BGR2RGB)
plt.imshow(flat_chess)
```




    <matplotlib.image.AxesImage at 0x7fd29ca4fdd0>




![png](output_3_1.png)



```python
# how to change gree nchessboard above to gray scale 
gray_flat_chess = cv2.cvtColor(flat_chess, cv2.COLOR_BGR2GRAY)
plt.imshow(gray_flat_chess, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7fd29b17ebd0>




![png](output_4_1.png)



```python
# how to upload anf get correct coloring for chessboard picture 
real_chess = cv2.imread("chessboard.jpg")
real_chess = cv2.cvtColor(real_chess, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7fd29b0fbe50>




![png](output_6_1.png)



```python
# how to turn chessboard picture into gray scale 
gray_real_chess = cv2.cvtColor(real_chess, cv2.COLOR_BGR2GRAY)
plt.imshow(gray_real_chess, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7fd29b07c250>




![png](output_7_1.png)



```python
# how to do Harris Corner Detection 
# block size is how far out it looks from a perspective corner 
# k size is aperature parameter 
gray = np.float32(gray_flat_chess)
dst = cv2.cornerHarris(src = gray, blockSize = 2, ksize = 3, k = 0.04)

dst = cv2. dilate(dst, None)
```


```python
# use our DST to see our flat chess board(green) make corners in it red 
flat_chess[dst>0.01*dst.max()]=[255,0,0]

plt.imshow(flat_chess)
```




    <matplotlib.image.AxesImage at 0x7fd29afe6fd0>




![png](output_9_1.png)



```python
# do what we did above but on the chessboard picture 
# catches where there are distinct changes in color that's why the shine on the pieces is red as well 
# colors the corners red bc of 255, can use green or blue as well instead of zero
gray = np.float32(gray_real_chess)
dst = cv2.cornerHarris(src = gray, blockSize = 2, ksize = 3, k = 0.04)
dst = cv2.dilate(dst, None)

real_chess[dst>0.01*dst.max()] = [255,0,0]

plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7fd29af5a910>




![png](output_10_1.png)



```python
# Shi-Tomasi Corner Detection 
# in parenthases: (which picture, max level of parameters we are detecting, quality level, minimum distance)
corners = cv2.goodFeaturesToTrack(gray_flat_chess, 64, 0.01, 10)
```


```python
corners = np.int0(corners)

for i in corners: 
    x,y = i.ravel()
    cv2.circle(flat_chess, (x,y), 3, (255,0,0,), -1)
    
plt.imshow(flat_chess)
```




    <matplotlib.image.AxesImage at 0x7fd29af54350>




![png](output_12_1.png)



```python
corners = cv2.goodFeaturesToTrack(gray_real_chess, 100, 0.01, 10)

corners = np.int0(corners)

for i in corners: 
    x,y = i.ravel()
    cv2.circle(real_chess, (x,y), 3, (0,255,0), -1)
    
plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7fd29aebc410>




![png](output_13_1.png)



```python

```
