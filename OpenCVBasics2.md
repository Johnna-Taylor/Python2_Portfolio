```python
# make sure to go in terminal and install: pip install --upgrade pip
# make sure to then go to terminal and install: pip install opencv-python
import cv2
import matplotlib.pyplot as plt 
%matplotlib inline
```


```python
# make sure to go to terminal and install: pip install 
import cv2
import matplotlib.pyplot as plt 
%matplotlib inline
```


```python
img = cv2.imread("mushroom.jpg")
```


```python
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7f207d7095d0>




![png](output_3_1.png)



```python
img1 = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7f207ce44510>




![png](output_5_1.png)



```python
# how to show image in hue saturation 
img2 = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
```


```python
plt.imshow(img2)
```




    <matplotlib.image.AxesImage at 0x7f207cdbc850>




![png](output_7_1.png)



```python
img3 = cv2.cvtColor(img, cv2.COLOR_BGR2HLS)
```


```python
plt.imshow(img3)
```




    <matplotlib.image.AxesImage at 0x7f207cd2b810>




![png](output_9_1.png)



```python
img1 = cv2.imread('do-not-copy-stamp.jpg')
img2 = cv2.imread('mushroom.jpg')
```


```python
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7f207cc97f50>




![png](output_11_1.png)



```python
img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)
img2 = cv2.cvtColor(img2, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7f207cc0db50>




![png](output_13_1.png)



```python
plt.imshow(img2)
```




    <matplotlib.image.AxesImage at 0x7f207cc0a250>




![png](output_14_1.png)



```python
img1 = cv2.resize(img1, (1200,1200))
img2 = cv2.resize(img2, (1200,1200))
```


```python
alpha = 0.5
beta = 0.5
```


```python
blended = cv2.addWeighted(img1, alpha, img2, beta, gamma=0)
```


```python
# creating this is really important for fluroescent microscopy 
plt.imshow(blended)
```




    <matplotlib.image.AxesImage at 0x7f207cb7a350>




![png](output_18_1.png)



```python
# how to make one seem way darker than the other 
alpha = 0.2
beta = 0.8

blended1 = cv2.addWeighted(img1, alpha, img2, beta, 0)
plt.imshow(blended1)
```




    <matplotlib.image.AxesImage at 0x7f207cad7b90>




![png](output_19_1.png)



```python
# this code and the one below it are how to put together/ in the corner of the other 
img1 = cv2.imread('do-not-copy-stamp.jpg')
img2 = cv2.imread('mushroom.jpg')

img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)
img2 = cv2.cvtColor(img2, cv2.COLOR_BGR2RGB)

img1 = cv2.resize(img1, (600,600))
```


```python
large_img = img2
small_img = img1

x_offset = 0 
y_offset = 0 

x_end = x_offset + small_img.shape[1]
y_end = y_offset + small_img.shape[0]

large_img[y_offset:y_end, x_offset:x_end] = small_img

plt.imshow(large_img)
```




    <matplotlib.image.AxesImage at 0x7f207cab7850>




![png](output_21_1.png)

