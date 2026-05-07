import numpy as np


TotalForce = np.array([0, 0, 0])
TotalTorque = np.array([0, 0, 0])


I0 = 1.5
w0 = 1.0


class RectangleLightSail():
    def __init__(self, width, height, resolution, StartingX, StartingY):
        
        self.width = width
        self.height = height
        self.resolution = resolution
        self.StartingX = StartingX
        self.StartingY = StartingY


    def compute(self):
        x = self.StartingX
        y = self.StartingY
        
        Force = np.array([0, 0, 0])
        Torque = np.array([0, 0, 0])
        
        dx = self.width / self.resolution
        dy = self.height / self.resolution
        
        dA = dx * dy
        for i in range(self.width):
            y = self.StartingY
            for j in range(self.height):
                r = np.sqrt(x**2 + y**2)
                I = I0 * np.exp(-2 * r**2 / w0**2)
                position = np.array([x, y, 0])
                
                normal = np.array([0, 0, 1])
                F = I * dA * normal
                
                tau = np.cross(position, F)
                
                Force = Force + F
                Torque = Torque + tau
                
                y += dy
            x += dx
        return Force, Torque

nyxRect1 = RectangleLightSail(10, 10, 50, -5, -5)
nyxRect2 = RectangleLightSail(10, 10, 50, -3, -5) #Different test case. The sail is placed assymtrically, relativaly to the centre is also close to the centre.

TotalForce, TotalTorque = nyxRect1.compute()

print(f"Test Case 1: \n ================== \n Total force is {TotalForce} \n Total torque is {TotalTorque} \n ==================")

TotalForce, TotalTorque = nyxRect2.compute()
print(f"Test Case 2: \n ================== \n Total force is {TotalForce} \n Total torque is {TotalTorque}")