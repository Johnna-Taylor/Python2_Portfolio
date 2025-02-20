# Python2_Portfolio
This is the portfolio of python that I learned during BISC450C-V84

## Feature Detection 
### Object Detection

In this class section we learned how to take an object from one picture and display how many of that object is in another picture with a heatmap 

```python
# make sure to go to terminal and install: python -m pip install -- upgrade pip
# then install: pip install opencv-python
import cv2
```


```python
import numpy as np 
```


```python
import matplotlib.pyplot as  plt
```


```python
%matplotlib inline
```


```python
# make sure to downlaod a photo of a singular sunflower and one of  field 
# uploading the single sunflower file 
full = cv2.imread('sunflower-training.jpg')
```


```python
# converting single sunflower image to color 
full = cv2.cvtColor(full, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(full)
```




    <matplotlib.image.AxesImage at 0x7fe59174ee10>



![output_6_1](https://github.com/user-attachments/assets/b3832e61-028b-4d87-bf22-27b722275bec)





```python
# upload the testing photo 
test = cv2.imread('sunflower-testing.jpg')
```


```python
test = cv2.cvtColor(test, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(test)
```




    <matplotlib.image.AxesImage at 0x7fe591683f10>




![output_9_1](https://github.com/user-attachments/assets/097d2ff8-1aed-4396-b76a-c560c038b409)




```python
# how to get the dimensions of both images 
print('Test image shape:', full.shape)
print('Training image shape:', test.shape)
```

    Test image shape: (600, 600, 3)
    Training image shape: (680, 1024, 3)



```python
methods = ['cv2.TM_CCOEFF', 'cv2.TM_CCOEFF_NORMED', 'cv2.TM_CCORR', 'cv2.TM_CCORR_NORMED', 'cv2.TM_SQDIFF', 'cv2.TM_SQDIFF_NORMED']
```


```python
# res is for template mapping
for m in methods: 
    
    test_copy = test.copy()
    method = eval(m)
    
    res = cv2.matchTemplate(test_copy,full, method)
    
    min_val, max_val, min_loc, max_loc = cv2.minMaxLoc(res)
    
# drawing the rectangle on our picture that we are matching 

    if method in[cv2.TM_SQDIFF, cv2.TM_SQDIFF_NORMED]:
        top_left = min_loc
    else: 
        top_left = max_loc
  
    height, width, channels = full.shape
    bottom_right = (top_left[0] + width, top_left[1] + height)
    
    cv2.rectangle(test_copy, top_left, bottom_right, (255,0,0), 10)
    
    # plotting 
    plt.subplot(121)
    plt.imshow(res)
    plt.title("Heatmap of template matching")
    plt.subplot(122)
    plt.imshow(test_copy)
    plt.title('Detection of template')

    #subtitle on a graph 
    plt.suptitle(m)
    
    # how to add a space/line between graphs that appear 
    plt.show()
    print('\n')
    print('\n')
```


![output_12_0](https://github.com/user-attachments/assets/ef2a769f-9046-4025-9798-3a93decd04c7)



    
    
    
    



![output_12_2](https://github.com/user-attachments/assets/6c5b2a03-20c4-4531-a185-700bb2f90275)



    
    
    
    



![output_12_4](https://github.com/user-attachments/assets/958ff37e-e516-4328-803c-514a396adf04)



    
    
    
    



![output_12_6](https://github.com/user-attachments/assets/31df0566-4b73-4cc6-aa66-75fe3a14ec61)



    
    
    
    



![output_12_8](https://github.com/user-attachments/assets/44f8d117-71a5-40ac-8df1-0fd12b068cf9)



    
    
    
    



![output_12_10](https://github.com/user-attachments/assets/dff2ad06-0deb-4cc3-aae6-5b456dab6572)


    
    
    
    



```python

```
