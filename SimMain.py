import numpy as np
import matplotlib.pyplot as plt


c = 3e8 #The Speed of Light in m/s. This is a constant for the simulation, and can be changed to test different scenarios.

I0 = 1e9 #The Intensity at The Exact Center of the Beam. This is a constant for the simulation, and can be changed to test different scenarios.
w0 = 1.0 #The beam waist, which is the radius at which the intensity falls to 1/e^2 of its maximum value. This is also a constant for the simulation, and can be changed to test different scenarios.

#Gaussian Beam Equation: I(r) = I0 * exp(-2 * r^2 / w0^2), How the intensity of the beam changes with distance from the center. This is a fundamental equation for simulating the interaction of light with the sail, and is used to calculate the force and torque on each element of the sail based on its position relative to the center of the beam.
def gaussian(r):
    return I0 * np.exp(-2 * r**2 / w0**2)

def rotateX(vector, angle):
    x, y, z = vector 
    
    return np.array([
        x,
        y*np.cos(angle) - z*np.sin(angle),
        y*np.sin(angle) + z*np.cos(angle)
    ])

def rotateY(vector, angle):
    x, y, z = vector
    
    return np.array([
        x*np.cos(angle) + z*np.sin(angle),
        y,
        -x*np.sin(angle) + z*np.cos(angle)
    ])

class RectangleLightSail():
    def __init__(self, width, height, resolution, StartingX=0, StartingY=0, ThetaX=0, ThetaY=0):
        self.width = width
        self.height = height
        self.resolution = resolution
        self.CenterX = StartingX
        self.CenterY = StartingY
        self.thetaX = ThetaX
        self.thetaY = ThetaY
        self.Force = None
        self.Torque = None
        self.xPoints=[]
        self.yPoints=[]
        self.zPoints=[]
        self.intensities = []
        self.forceVectors = []

    def compute(self):
        self.xPoints=[]
        self.yPoints=[]
        self.zPoints=[]
        self.intensities = []
        
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
                
                position = np.array([x, y, 0.0])
                
                position = rotateX(position, self.thetaX)
                position = rotateY(position, self.thetaY)
                
                r = np.sqrt(position[0]**2 + position[1]**2)
                I = gaussian(r)
                
                self.xPoints.append(position[0])
                self.yPoints.append(position[1])
                self.zPoints.append(position[2])
                self.intensities.append(I)
                
                normal = np.array([np.sin(self.thetaY) * np.cos(self.thetaX), -np.sin(self.thetaX), np.cos(self.thetaX) * np.cos(self.thetaY)])
                
                cosTheta = np.dot(normal, np.array([0, 0, 1]))
                
                if cosTheta > 0:
                    fmag = (2 *I / 3e8) * (cosTheta**2) * dA
                    
                    F = fmag * normal
                else:
                    F = np.array([0.0, 0.0, 0.0])
                
                self.forceVectors.append(F)
                
                tau = np.cross(position, F)
                
                Force = Force + F
                Torque = Torque + tau
        self.Force = Force
        self.Torque = Torque
        return Force, Torque

    def Visualize(self):
        fig = plt.figure()
        
        Force, Torque = self.compute()
        
        ax = fig.add_subplot(111, projection='3d')
        scatter = ax.scatter(self.xPoints, self.yPoints, self.zPoints, c=self.intensities, cmap='plasma', s=30)
        
        fig.colorbar(scatter, ax=ax, label='Beam Intensity (W/m^2)')
        
        ax.quiver(0, 0, 0, Force[0], Force[1], Force[2], color='red', length=2, normalize=True)
        
        ax.text2D(
            0.05, 0.95, f'Total Force: {Force}\nTotal Torque: {Torque}'
        )
        
        ax.set_box_aspect([1,1,1])
        ax.view_init(elev=30, azim=45)
        
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        
        # ax.margins(0.5)
        
        plt.show()

class SphereLightSail():
    def __init__(self, radius, resolution, thetaX, thetaY):
        self.radius = radius
        self.resolution = resolution
        self.thetaX = thetaX
        self.thetaY = thetaY
        self.Force = None
        self.Torque = None
        self.xPoints=[]
        self.yPoints=[]
        self.zPoints=[]
        self.intensities = []

    def compute(self):
        self.xPoints=[]
        self.yPoints=[]
        self.zPoints=[]
        
        Force = np.array([0.0, 0.0, 0.0])
        Torque = np.array([0.0, 0.0, 0.0])
        
        for i in range(self.resolution):
            for j in range(self.resolution):
                
                theta = (np.pi / 2) * (i + 0.5) / self.resolution
                phi = 2 * np.pi * (j + 0.5) / self.resolution
                
                x = self.radius * np.sin(theta) * np.cos(phi)
                y = self.radius * np.sin(theta) * np.sin(phi)
                z = self.radius * np.cos(theta)
                
                position = np.array([x, y, z])
                
                position = rotateX(position, self.thetaX)
                position = rotateY(position, self.thetaY)
                
                self.xPoints.append(position[0])
                self.yPoints.append(position[1])
                self.zPoints.append(position[2])
                
                normal = position / np.linalg.norm(position)
                
                r = np.sqrt(position[0]**2 + position[1]**2)
                I = gaussian(r)
                
                self.intensities.append(I)
                
                dtheta = (np.pi / 2) / self.resolution
                dphi = (2 * np.pi) / self.resolution
                
                dA = self.radius**2 * np.sin(theta) * dtheta * dphi
                
                # F = I * dA * normal
                
                cosTheta = np.dot(normal, np.array([0, 0, 1]))
                
                if cosTheta > 0:
                    fmag = (2 *I / c) * (cosTheta**2) * dA
                    
                    F = fmag * normal
                else:
                    continue
                
                tau = np.cross(position, F)
                
                Force = Force + F
                Torque = Torque + tau
        self.Force = Force
        self.Torque = Torque
        return Force, Torque
        
    def Visualize(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        scatter = ax.scatter(self.xPoints, self.yPoints, self.zPoints, c=self.intensities, cmap='plasma', s=30)
        
        fig.colorbar(scatter, ax=ax, label='Beam Intensity (W/m^2)')
        
        ax.quiver(0, 0, 0, self.Force[0], self.Force[1], self.Force[2], color='red', length=2, normalize=True)
        
        ax.set_box_aspect([1,1,1])
        ax.view_init(elev=30, azim=45)
        
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        
        ax.margins(0.5)
        
        plt.show()

class ParaboloidLightSail(): #Paraboloid Reflector
    def __init__(self, radius, a, resolution, thetaX, thetaY):
        self.radius = radius
        self.a = a
        self.resolution = resolution
        self.thetaX = thetaX
        self.thetaY = thetaY
        self.Force = None
        self.Torque = None
        self.xPoints=[]
        self.yPoints=[]
        self.zPoints=[]
        self.intensities = []
    
    def compute(self):
        self.xPoints=[]
        self.yPoints=[]
        self.zPoints=[]
        self.intensities = []
        
        Force = np.array([0.0, 0.0, 0.0])
        Torque = np.array([0.0, 0.0, 0.0])
        
        for i in range(self.resolution):
            for j in range(self.resolution):
                r = self.radius * (i + 0.5) / self.resolution
                phi = 2 * np.pi * (j + 0.5) / self.resolution
                
                x = r * np.cos(phi)
                y = r * np.sin(phi)
                
                z = -self.a * (self.radius**2 - r**2)
                
                normal = np.array([2 * self.a * x, 2 * self.a * y, 1.0])
                normal = normal / np.linalg.norm(normal)
                
                normal = rotateX(normal, self.thetaX)
                normal = rotateY(normal, self.thetaY)
                
                position = np.array([x, y, z])
                
                position = rotateX(position, self.thetaX)
                position = rotateY(position, self.thetaY)
                
                self.xPoints.append(position[0])
                self.yPoints.append(position[1])
                self.zPoints.append(position[2])
                
                beam = np.sqrt(position[0]**2+position[1]**2)
                
                I = gaussian(beam)
                
                dr = self.radius / self.resolution
                dphi = 2 * np.pi / self.resolution
                
                dA = r * dr * dphi
                
                cosTheta = np.dot(normal, np.array([0, 0, 1]))
                
                if cosTheta > 0:
                    fmag = (2 *I / c) * (cosTheta**2) * dA
                    
                    F = fmag * normal
                else:
                    continue
                
                tau = np.cross(position, F)
                
                Force =+ F
                Torque =+ tau
                
                self.intensities.append(I)
        self.Force=Force
        self.Torque=Torque
        return Force, Torque
        
    
    def Visualize(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        scatter = ax.scatter(self.xPoints, self.yPoints, self.zPoints, c=self.intensities, cmap='plasma', s=30)
        
        fig.colorbar(scatter, ax=ax, label='Beam Intensity (W/m^2)')
        
        ax.quiver(0, 0, 0, self.Force[0], self.Force[1], self.Force[2], color='red', length=2, normalize=True)
        
        ax.set_box_aspect([1,1,1])
        ax.view_init(elev=30, azim=45)
        
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        
        ax.margins(0.5)
        
        plt.show()

nyxRect1 = RectangleLightSail(10, 10, 50, 0, 0, np.radians(0), np.radians(0))
nyxRect2 = RectangleLightSail(10, 10, 50, 0, 0, np.radians(20), np.radians(0))
nyxRect3 = RectangleLightSail(10, 10, 50, 0, 0, np.radians(45), np.radians(0))
nyxRect4 = RectangleLightSail(10, 10, 50, 0, 0, np.radians(90), np.radians(0))

nyxSphere1 = SphereLightSail(1, 100, np.radians(0), np.radians(0))
nyxSphere2 = SphereLightSail(1, 100, np.radians(30), np.radians(0))
nyxSphere3 = SphereLightSail(1, 100, np.radians(60), np.radians(0))

nyxParaboloid1 = ParaboloidLightSail(1, 0.5, 100, 0, 0)

TotalForce, TotalTorque = nyxParaboloid1.compute()
print(f"Paraboloid Reflector Test Case 1: \n ================== \n Total force is {TotalForce} \n Total torque is {TotalTorque} \n ==================")
nyxParaboloid1.Visualize()

# TotalForce, TotalTorque = nyxSphere1.compute()
# print(f"Sphere Test Case 1: \n ================== \n Total force is {TotalForce} \n Total torque is {TotalTorque} \n ==================")
# nyxSphere1.Visualize()

# TotalForce, TotalTorque = nyxSphere2.compute()
# print(f"Sphere Test Case 2: \n ================== \n Total force is {TotalForce} \n Total torque is {TotalTorque} \n ==================")
# nyxSphere2.Visualize()

# TotalForce, TotalTorque = nyxSphere3.compute()
# print(f"Sphere Test Case 3: \n ================== \n Total force is {TotalForce} \n Total torque is {TotalTorque} \n ==================")
# nyxSphere3.Visualize()


# TotalForce, TotalTorque = nyxRect1.compute()
# print(f"Test Case 1: \n ================== \n Total force is {TotalForce} \n Total torque is {TotalTorque} \n ==================")
# nyxRect1.Visualize()

# TotalForce, TotalTorque = nyxRect2.compute()
# print(f"Test Case 2: \n ================== \n Total force is {TotalForce} \n Total torque is {TotalTorque} \n ==================")

# TotalForce, TotalTorque = nyxRect3.compute()
# print(f"Test Case 3: \n ================== \n Total force is {TotalForce} \n Total torque is {TotalTorque} \n ==================")

# TotalForce, TotalTorque = nyxRect4.compute()
# print(f"Test Case 4: \n ================== \n Total force is {TotalForce} \n Total torque is {TotalTorque} \n ==================")