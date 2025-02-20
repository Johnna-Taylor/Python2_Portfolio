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




![png](output_6_1.png)



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




![png](output_9_1.png)



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


![png](output_12_0.png)


    
    
    
    



![png](output_12_2.png)


    
    
    
    



![png](output_12_4.png)


    
    
    
    



![png](output_12_6.png)


    
    
    
    



![png](output_12_8.png)


    
    
    
    



![png](output_12_10.png)


    
    
    
    



```python

```
