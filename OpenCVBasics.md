# Python2_Portfolio
This is the portfolio of python that I learned during BISC450C-V84

## Open CV
### Part 1

In this class section we learned how to import a photo, display it in the correct colors, scale the photo, and rotating the photo

```python
# make sure to go to terminal and do: pip install --upgrade pip 
# then do: pip install opencv-python
# and make sure to download a picture: example is mushroom but can do dog (jpg file only)
import numpy as np 
import matplotlib.pyplot as plt 
%matplotlib inline
```


```python
import cv2
```


```python
img = cv2.imread("mushroom.jpg")
```


```python
type(img)
```




    numpy.ndarray




```python
# can use this to make sure our image was imported correctly 
img_wrong = cv2.imread('wrong/path/doesnot/abcdegh.jpg')
```


```python
type(img_wrong)
```




    NoneType




```python
# matplotlib expects a different order of colors, that's why the mushrooms are blue and not red like the original picture 
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7f4584081110>




![output_6_1](https://github.com/user-attachments/assets/df218f21-1ec8-4551-bd8f-2b40babc660b)




```python
# how to switch it from reading colors blue/green/red to red/green/bue 
fix_img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(fix_img)
```




    <matplotlib.image.AxesImage at 0x7f4566fb2550>




![output_8_1](https://github.com/user-attachments/assets/aa9d1535-57ea-41e6-a4f9-a593e0c7d60c)




```python
# can show the image on a gray scale 
# shows the dimensions of the picture 
img_gray = cv2.imread("mushroom.jpg", cv2.IMREAD_GRAYSCALE)
img_gray.shape
```




    (1929, 3000)




```python
plt.imshow(img_gray)
```




    <matplotlib.image.AxesImage at 0x7f4566f320d0>




![output_10_1](https://github.com/user-attachments/assets/88586db3-06d3-4a4f-92b9-89070931746b)




```python
# those colors because of matplotlib so we must transform it again
plt.imshow(img_gray, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7f4566e9e2d0>




![output_11_1](https://github.com/user-attachments/assets/34330002-5613-40de-a6b2-7a7a3a209ea4)




```python
# we can resize our images 
fix_img.shape
```




    (1929, 3000, 3)




```python
new_img = cv2.resize(fix_img,(2300,2500))
plt.imshow(new_img)
```




    <matplotlib.image.AxesImage at 0x7f4566e0b990>




![output_13_1](https://github.com/user-attachments/assets/154f1ae4-6647-468a-8e13-916c5ee6b7c9)




```python
new_img.shape
```




    (2500, 2300, 3)




```python
# can use a ratio to scale the picture 
w_ratio = 0.5 
h_ratio = 0.5 

new_img = cv2.resize(fix_img, (0,0), fix_img, w_ratio, h_ratio)
```


```python
plt.imshow(new_img)
```




    <matplotlib.image.AxesImage at 0x7f4566dedb10>




![output_16_1](https://github.com/user-attachments/assets/9b69e9f8-16ec-4211-be61-c049ee9d6d3f)




```python
new_img.shape
```




    (964, 1500, 3)




```python
# we can flip images: 0 just means upside down: -1 means upside down and horizontally flipped 
flip_img = cv2.flip(fix_img, 0)
plt.imshow(flip_img)
```




    <matplotlib.image.AxesImage at 0x7f4566d5ad50>




![output_18_1](https://github.com/user-attachments/assets/3c992246-e4eb-4472-a97a-f8c8e2e58e71)




```python
flip_img2 = cv2.flip(fix_img, -1)
plt.imshow(flip_img2)
```




    <matplotlib.image.AxesImage at 0x7f4566cc9e90>




![output_19_1](https://github.com/user-attachments/assets/3a46b12e-d1a5-4c0d-b969-0e5be50afe09)




```python
type(fix_img)
```




    numpy.ndarray




```python
# how to save an image: only saves the flip NOT the color conversion 
cv2.imwrite('mushroom_fixed_image.jpg', flip_img)
```




    True


