import numpy as np
import matplotlib.pyplot as plt


I0 = 1.5 #The Intensity at The Exact Center of the Beam. This is a constant for the simulation, and can be changed to test different scenarios.
w0 = 1.0 #The beam waist, which is the radius at which the intensity falls to 1/e^2 of its maximum value. This is also a constant for the simulation, and can be changed to test different scenarios.

#Gaussian Beam Equation: I(r) = I0 * exp(-2 * r^2 / w0^2), How the intensity of the beam changes with distance from the center. This is a fundamental equation for simulating the interaction of light with the sail, and is used to calculate the force and torque on each element of the sail based on its position relative to the center of the beam.
def gaussian(r):
    return I0 * np.exp(-2 * r**2 / w0**2)

class RectangleLightSail():
    def __init__(self, width, height, resolution, StartingX=0, StartingY=0, ThetaX=0, ThetaY=0):
        self.width = width
        self.height = height
        self.resolution = resolution
        self.CenterX = StartingX
        self.CenterY = StartingY
        self.thetaX = ThetaX
        self.thetaY = ThetaY

    def compute(self):
        x = self.CenterX
        y = self.CenterY
        
        Force = np.array([0.0, 0.0, 0.0])
        Torque = np.array([0.0, 0.0, 0.0])
        
        dx = self.width / self.resolution
        dy = self.height / self.resolution
        
        dA = dx * dy
        for i in range(self.resolution):
            for j in range(self.resolution):
                
                x = (self.CenterX - self.width / 2) + (i+0.5) * dx
                y = (self.CenterY - self.height / 2) + (j+0.5) * dy
                
                r = np.sqrt(x**2 + y**2)
                I = gaussian(r)
                position = np.array([x, y, 0.0])
                
                normal = np.array([np.sin(self.thetaY) * np.cos(self.thetaX), -np.sin(self.thetaX), np.cos(self.thetaX) * np.cos(self.thetaY)])
                # print(f"Theta_X: {self.thetaX}, Normal Vector: {normal}")
                F = I * dA * normal
                
                tau = np.cross(position, F)
                
                Force = Force + F
                Torque = Torque + tau
        # print(normal)
        return Force, Torque

class SphereLightSail():
    def __init__(self, radius, resolution):
        self.radius = radius
        self.resolution = resolution
        self.xPoints=[]
        self.yPoints=[]
        self.zPoints=[]

    def generatePatches(self):
        self.xPoints=[]
        self.yPoints=[]
        self.zPoints=[]
        for i in range(self.resolution):
            for j in range(self.resolution):
                
                theta = (np.pi / 2) * (i + 0.5) / self.resolution
                phi = 2 * np.pi * (j + 0.5) / self.resolution
                
                x = self.radius * np.sin(theta) * np.cos(phi)
                y = self.radius * np.sin(theta) * np.sin(phi)
                z = self.radius * np.cos(theta)
                
                position = np.array([x, y, z])
                
                self.xPoints.append(x)
                self.yPoints.append(y)
                self.zPoints.append(z)
                
                normal = position / self.radius
            # print(f"Patch {i*self.resolution + j}: \nPosition: {position}, \nNormal: {normal} \n")
    def VerifySphere(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(self.xPoints, self.yPoints, self.zPoints)
        
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        
        plt.show()


nyxRect1 = RectangleLightSail(10, 10, 50, 0, 0, np.radians(0), np.radians(0))
nyxRect2 = RectangleLightSail(10, 10, 50, 0, 0, np.radians(0), np.radians(20))
nyxRect3 = RectangleLightSail(10, 10, 50, 0, 0, np.radians(45), np.radians(0))
nyxRect4 = RectangleLightSail(10, 10, 50, 0, 0, np.radians(90), np.radians(90))

nyxSphere1 = SphereLightSail(1, 50)
nyxSphere1.generatePatches()
nyxSphere1.VerifySphere()

# TotalForce, TotalTorque = nyxRect1.compute()
# print(f"Test Case 1: \n ================== \n Total force is {TotalForce} \n Total torque is {TotalTorque} \n ==================")

# TotalForce, TotalTorque = nyxRect2.compute()
# print(f"Test Case 2: \n ================== \n Total force is {TotalForce} \n Total torque is {TotalTorque} \n ==================")

# TotalForce, TotalTorque = nyxRect3.compute()
# print(f"Test Case 3: \n ================== \n Total force is {TotalForce} \n Total torque is {TotalTorque} \n ==================")

# TotalForce, TotalTorque = nyxRect4.compute()
# print(f"Test Case 4: \n ================== \n Total force is {TotalForce} \n Total torque is {TotalTorque} \n ==================")