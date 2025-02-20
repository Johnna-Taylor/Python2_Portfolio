```python
# make sure to upgrade pip and pip install opencv-python
import cv2
import numpy as np
import matplotlib.pyplot as plt 
%matplotlib inline
```


```python
# download a picture of lots of cereal box logos and then another picture with just one of the matching box logos 
def display(img, cmap = 'gray'):
    fig = plt.figure(figsize = (12,10))
    ax = fig.add_subplot(111)
    ax.imshow(img, cmap = 'gray')
```


```python
apple_jacks = cv2.imread("apple-jacks.jpg", 0)
display(apple_jacks)
```


![png](output_2_0.png)



```python
cereals = cv2.imread("many-cereals.jpg", 0)
display(cereals)
```


![png](output_3_0.png)



```python
# a way of finding key points and descripters in the thing that you are searching
orb = cv2.ORB_create()

kp1,des1 = orb.detectAndCompute(apple_jacks, mask=None)
kp2,des2 = orb.detectAndCompute(cereals, mask=None)
```


```python
# were the matching actually takes place (bf: brute force) 
bf = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck = True)
matches = bf.match(des1, des2)
```


```python
# there will be matches that are not as good, so we write this line so that we only get the good matches 
matches = sorted(matches, key = lambda x:x.distance)
```


```python
apple_jacks_matches = cv2.drawMatches(apple_jacks, kp1, cereals, kp2, matches[:25], None, flags = 2)
```


```python
display(apple_jacks_matches)
```


![png](output_8_0.png)



```python
# go to terminal and do: pip uninstall opencv-python
# then do: pip install oprncv-contrib-python
sift = cv2.SIFT_create()
```


```python
kp1, des1 = sift.detectAndCompute(apple_jacks, None)
kp2, des2 = sift.detectAndCompute(cereals, None)
```


```python
# knn = no nearest neighbor
# it finds the best matches for the descriptors for this set 
bf = cv2.BFMatcher()
matches = bf.knnMatch(des1, des2, k=2)
```


```python
# applying a ratio test 
# if the matches are distanced apart thn it's not a good match, if the matches are clustered together then we probably have the correct box 
good = []

for match1, match2 in matches: 
    if match1.distance < 0.75*match2.distance:
        good.append([match1])
```


```python
print('Length of total matches:', len(matches))
print('Length of good matches:', len(good))
```

    Length of total matches: 5134
    Length of good matches: 99



```python
sift_matches = cv2.drawMatchesKnn(apple_jacks, kp1, cereals, kp2, good, None, flags =2)
display(sift_matches)
```


![png](output_14_0.png)



```python
# flan based matching 
# it doen't find the best mathes but general good matches 
sift = cv2.SIFT_create()

kp1, des1 = sift.detectAndCompute(apple_jacks, None)
kp2, des2 = sift.detectAndCompute(cereals, None)
```


```python
# defining parameters 
flann_index_KDtree = 0 
index_params = dict(algorithm=flann_index_KDtree, trees=5)
search_params = dict(checks=50)
```


```python
# running the matchers 
flann = cv2.FlannBasedMatcher(index_params, search_params)

matches = flann.knnMatch(des1, des2, k=2)

good = []

for match1, match2, in matches: 
    if match1.distance < 0.75*match2.distance: 
        good.append([match1])
```


```python
# draw lines of matches from one picture to the other
flann_matches = cv2.drawMatchesKnn(apple_jacks, kp1, cereals, kp2, good, None, flags = 0)
display(flann_matches)
```


![png](output_18_0.png)



```python
# adding a mask to it 
sift = cv2.SIFT_create()

kp1, des1 = sift.detectAndCompute(apple_jacks, None)
kp2, des2 = sift.detectAndCompute(cereals, None)
```


```python
flann_index_KDtree = 0
index_params = dict(algorithm=flann_index_KDtree, trees = 5)
search_params = dict(check=50)
```


```python
flann = cv2.FlannBasedMatcher(index_params, search_params)

matches = flann.knnMatch(des1, des2, k=2)
```


```python
# 0:0 means that pure black will be added 
matchesMask = [[0,0] for i in range(len(matches))]
```


```python
# how to display certain points on a cereal box and use another color to connect the two to see where it came from 
for i, (match1, match2) in enumerate(matches):
    if match1.distance < 0.75*match2.distance:
        matchesMask[i]=[1,0]
        
draw_params = dict(matchColor = (0,255,0),
                  singlePointColor = (255,0,0),
                  matchesMask = matchesMask, 
                  flags =0)
```


```python
flann_matches = cv2.drawMatchesKnn(apple_jacks, kp1, cereals, kp2, matches, None, **draw_params)
display(flann_matches)
```


![png](output_24_0.png)

