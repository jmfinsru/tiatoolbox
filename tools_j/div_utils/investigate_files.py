import numpy as np

arr = np.load("/home/jenny/pannuke_dataset/Fold_3/masks/masks.npy")

# Basic info
print("Type:", type(arr))
print("Shape:", arr.shape)
print("Data type:", arr.dtype)

# First 3 rows/elements
print(arr[:3])  
