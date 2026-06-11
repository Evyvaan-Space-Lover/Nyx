import numpy as np


TotalForce = np.array([0, 0, 0])
TotalTorque = np.array([0, 0, 0])


I0 = 1.5 #The Intensity at The Exact Center of the Beam. This is a constant for the simulation, and can be changed to test different scenarios.
w0 = 1.0 #The beam waist, which is the radius at which the intensity falls to 1/e^2 of its maximum value. This is also a constant for the simulation, and can be changed to test different scenarios.

#Gaussian Beam Equation: I(r) = I0 * exp(-2 * r^2 / w0^2), How the intensity of the beam changes with distance from the center. This is a fundamental equation for simulating the interaction of light with the sail, and is used to calculate the force and torque on each element of the sail based on its position relative to the center of the beam.
def gaussian(r):
    return I0 * np.exp(-2 * r**2 / w0**2)

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
        for i in range(self.resolution):
            for j in range(self.resolution):
                
                x = (self.StartingX - self.width / 2) + (i+0.5) * dx
                y = (self.StartingY - self.height / 2) + (j+0.5) * dy

                r = np.sqrt(x**2 + y**2)
                I = gaussian(r)
                position = np.array([x, y, 0])
                
                normal = np.array([0, 0, 1])
                F = I * dA * normal
                
                tau = np.cross(position, F)
                
                Force = Force + F
                Torque = Torque + tau
                
                
            
        return Force, Torque

nyxRect1 = RectangleLightSail(10, 10, 50, 0, 0)
nyxRect2 = RectangleLightSail(10, 10, 100, 0, 0)
nyxRect3 = RectangleLightSail(10, 10, 500, 0, 0)

TotalForce, TotalTorque = nyxRect1.compute()

print(f"Test Case 1: \n ================== \n Total force is {TotalForce} \n Total torque is {TotalTorque} \n ==================")

TotalForce, TotalTorque = nyxRect2.compute()
print(f"Test Case 2: \n ================== \n Total force is {TotalForce} \n Total torque is {TotalTorque}")