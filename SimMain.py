import numpy as np


TotalForce = np.array([0, 0, 0])
TotalTorque = np.array([0, 0, 0])


I0 = 1.5
w0 = 1.0


class RectangleLightSail():
    def __init__(self, width, height, StartingX, StartingY):
        
        self.width = width * 100
        self.height = height * 100
        self.StartingX = StartingX
        self.StartingY = StartingY


    def compute(self):
        x = self.StartingX
        y = self.StartingY
        
        Force = np.array([0, 0, 0])
        Torque = np.array([0, 0, 0])
        
        for i in range(self.width):
            y = self.StartingY
            for j in range(self.height):
                r = np.sqrt(x**2 + y**2)
                I = I0 * np.exp(-2 * r**2 / w0**2)
                position = np.array([x, y, 0])
                
                normal = np.array([0, 0, 1])
                F = I * normal
                
                tau = np.cross(position, F)
                
                Force = Force + F
                Torque = Torque + tau
                
                y += 1
            x += 1
        return Force, Torque

nyx1 = RectangleLightSail(10, 10, -5, -5)
TotalForce, TotalTorque = nyx1.compute()

print(f"total force is {TotalForce} \n total torque is {TotalTorque}")
