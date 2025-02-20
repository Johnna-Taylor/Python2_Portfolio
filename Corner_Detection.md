# Python2_Portfolio
This is the portfolio of python that I learned during BISC450C-V84

## Aspect Detection 
### Corner Detection

In this class section we learned how to put markers(colored dots) on corners in a photo

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




![output_3_1](https://github.com/user-attachments/assets/e5400137-a693-476f-bc2b-d8b9288b48d6)




```python
# how to change gree nchessboard above to gray scale 
gray_flat_chess = cv2.cvtColor(flat_chess, cv2.COLOR_BGR2GRAY)
plt.imshow(gray_flat_chess, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7fd29b17ebd0>



![output_4_1](https://github.com/user-attachments/assets/202fd0c7-db78-449d-a7b6-8863b24b8704)





```python
# how to upload anf get correct coloring for chessboard picture 
real_chess = cv2.imread("chessboard.jpg")
real_chess = cv2.cvtColor(real_chess, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7fd29b0fbe50>




![output_6_1](https://github.com/user-attachments/assets/15250f54-1268-4c0d-b531-98e051d107ae)




```python
# how to turn chessboard picture into gray scale 
gray_real_chess = cv2.cvtColor(real_chess, cv2.COLOR_BGR2GRAY)
plt.imshow(gray_real_chess, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7fd29b07c250>




![output_7_1](https://github.com/user-attachments/assets/827b8a9d-7bc6-47bf-b481-73cd1e0f2896)




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




![output_9_1](https://github.com/user-attachments/assets/0d2ffa6c-9ffb-4c69-a3d6-6c4691ec1b03)




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




![output_10_1](https://github.com/user-attachments/assets/e6d8db32-6421-403e-ba61-8d261333953b)




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




![output_12_1](https://github.com/user-attachments/assets/b95922ee-e47b-4e17-895f-daf82e6de675)




```python
corners = cv2.goodFeaturesToTrack(gray_real_chess, 100, 0.01, 10)

corners = np.int0(corners)

for i in corners: 
    x,y = i.ravel()
    cv2.circle(real_chess, (x,y), 3, (0,255,0), -1)
    
plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7fd29aebc410>




![output_13_1](https://github.com/user-attachments/assets/aef4854d-7e57-4ecc-9b13-b13815314fd7)




```python

```
