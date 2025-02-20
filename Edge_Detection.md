# Python2_Portfolio
This is the portfolio of python that I learned during BISC450C-V84

## Aspect Detection 
### Edge Detection

In this class section we learned how to line the edges of an image with dots and using a blurred threshold 

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




![output_2_1](https://github.com/user-attachments/assets/ac65bc27-d9ba-4e01-871c-22d89da20ec0)




```python
# setting the color ratio for the output pictures; 25 is the max but we want two colors 
edges = cv2.Canny(image =img, threshold1 = 127, threshold2 = 127)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7f2652083d90>




![output_3_1](https://github.com/user-attachments/assets/8e57b51e-24bb-4af3-96f0-566dc11a63e9)




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




![output_5_1](https://github.com/user-attachments/assets/08b9d8c2-2df3-4176-ace3-dd6f853ea9a3)




```python
edges = cv2.Canny(image = img, threshold1 = lower, threshold2 = upper + 100)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7f2651f6fd50>




![output_6_1](https://github.com/user-attachments/assets/d98517d1-7660-42a2-9882-a9369809ee4d)




```python
# blurring as a solution; creates lower resolution with pixels  
blurred_img = cv2.blur(img, ksize = (5,5))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower, 
                 threshold2 = upper)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7f2651ee2550>




![output_7_1](https://github.com/user-attachments/assets/396345b9-2d12-4ee0-9e94-9931e513ae0e)




```python
blurred_img = cv2.blur(img, ksize = (7,7))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower, 
                 threshold2 = upper)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7f2651ec98d0>




![output_8_1](https://github.com/user-attachments/assets/151b04b3-2b46-44c1-8126-a8972c0133f2)




```python
# can increase the upper threshold limit 
blurred_img = cv2.blur(img, ksize = (5,5))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower, 
                 threshold2 = upper + 50)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7f2651e3a990>




![output_9_1](https://github.com/user-attachments/assets/ddff902b-a95f-4e30-bcef-dca1d697b946)




```python
blurred_img = cv2.blur(img, ksize = (5,5))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower, 
                 threshold2 = upper + 100)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7f2651da7a90>




![output_10_1](https://github.com/user-attachments/assets/5ac17521-b885-47b4-ad6f-391a5f589ac8)




```python
blurred_img = cv2.blur(img, ksize = (7,7))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower, 
                 threshold2 = upper + 60)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7f2651d16bd0>




![output_11_1](https://github.com/user-attachments/assets/3f97b939-0687-4511-89e9-982236ec1e98)




```python
blurred_img = cv2.blur(img, ksize = (7,7))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower, 
                 threshold2 = upper + 60)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7f2651d03b50>




![output_12_1](https://github.com/user-attachments/assets/df7919b4-9f76-4e4f-b632-5f65446a5333)




```python
# I was just messing around trying to see how well i could get mine to outline
blurred_img = cv2.blur(img, ksize = (1,1))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower, 
                 threshold2 = upper + 145)
                  
plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7f2651c6bc90>




![output_13_1](https://github.com/user-attachments/assets/c318bc73-5c79-452b-b522-2ac7561c89dc)


