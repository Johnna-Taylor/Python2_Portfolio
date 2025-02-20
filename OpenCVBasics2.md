# Python2_Portfolio
This is the portfolio of python that I learned during BISC450C-V84

## Open CV 
### Part 2

In this class section we learned how to do hue saturation on a photo, as well as how to overlay one photo with another 

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




![output_3_1](https://github.com/user-attachments/assets/1560fa5f-7011-4858-aeeb-4d157d137b6f)




```python
img1 = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7f207ce44510>




![output_5_1](https://github.com/user-attachments/assets/f9b42a96-1ac6-4545-a698-abb0c1ad4f60)




```python
# how to show image in hue saturation 
img2 = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
```


```python
plt.imshow(img2)
```




    <matplotlib.image.AxesImage at 0x7f207cdbc850>




![output_7_1](https://github.com/user-attachments/assets/f9e5d69a-a7e8-47d0-bc52-b8e365d732be)




```python
img3 = cv2.cvtColor(img, cv2.COLOR_BGR2HLS)
```


```python
plt.imshow(img3)
```




    <matplotlib.image.AxesImage at 0x7f207cd2b810>



![output_9_1](https://github.com/user-attachments/assets/5c66991b-2753-4263-a5e4-4eeeeb539c0f)





```python
img1 = cv2.imread('do-not-copy-stamp.jpg')
img2 = cv2.imread('mushroom.jpg')
```


```python
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7f207cc97f50>




![output_11_1](https://github.com/user-attachments/assets/fab61983-8fce-40ad-a071-1be96ea3baa7)




```python
img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)
img2 = cv2.cvtColor(img2, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7f207cc0db50>




![output_13_1](https://github.com/user-attachments/assets/efbbbcca-5455-481b-a8b6-d571367654b0)




```python
plt.imshow(img2)
```




    <matplotlib.image.AxesImage at 0x7f207cc0a250>




![output_14_1](https://github.com/user-attachments/assets/13e77814-2c05-4f5e-863c-e8eb383269d0)




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




![output_18_1](https://github.com/user-attachments/assets/f0b47c69-1097-4680-9773-6a8dbcb7c9b9)




```python
# how to make one seem way darker than the other 
alpha = 0.2
beta = 0.8

blended1 = cv2.addWeighted(img1, alpha, img2, beta, 0)
plt.imshow(blended1)
```




    <matplotlib.image.AxesImage at 0x7f207cad7b90>




![output_19_1](https://github.com/user-attachments/assets/4087e603-7e30-4b4e-b818-7137f34489ef)




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




![output_21_1](https://github.com/user-attachments/assets/6ff11027-1d8b-40fa-affd-1812f5452114)


