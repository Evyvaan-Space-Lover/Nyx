import numpy as np

x = 0.5
y = 0.5
z = 0

I0 = 1.5
w0 = 1.0

r = np.sqrt(x**2 + y**2)

I = I0 * np.exp(-2 * r**2 / w0**2)

position = np.array([x, y, z])
print("Position: ", position)
print("Intensity: ", I)