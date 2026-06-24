import numpy as np
import matplotlib.pyplot as plt


c = 3e8 #The Speed of Light in m/s. This is a constant for the simulation, and can be changed to test different scenarios.

I0 = 1e9 #The Intensity at The Exact Center of the Beam. This is a constant for the simulation, and can be changed to test different scenarios.
w0 = 1.0 #The beam waist, which is the radius at which the intensity falls to 1/e^2 of its maximum value. This is also a constant for the simulation, and can be changed to test different scenarios.

beamX = 0
beamY = 0
beamDih = np.array([0.0, 0.0, 1.0]) #Beam direction, i thought it would be funny to name it "Beam Dih"
beamDih = beamDih / np.linalg.norm(beamDih)

#Gaussian Beam Equation: I(r) = I0 * exp(-2 * r^2 / w0^2), How the intensity of the beam changes with distance from the center. This is a fundamental equation for simulating the interaction of light with the sail, and is used to calculate the force and torque on each element of the sail based on its position relative to the center of the beam.
def gaussian(x, y):
    r = np.sqrt((x-beamX)**2 + (y-beamY)**2)
    
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
def calcForce(I, normal, dA):
    cosTheta = np.dot(normal, beamDih)
    if cosTheta > 1e-10:
        fmag = (I / c) * cosTheta * dA
        F = -fmag * (beamDih - 2*cosTheta*normal)
    else:
        F = np.array([0.0, 0.0, 0.0])
    
    return F

class RectangleLightSail():
    def __init__(self, width, height, resolution, StartingX=0, StartingY=0, ThetaX=0, ThetaY=0):
        self.width = width
        self.height = height
        self.resolution = resolution
        self.CenterX = StartingX
        self.CenterY = StartingY
        self.thetaX = ThetaX
        self.thetaY = ThetaY
        
        self.mass = 1.0
        
        self.velocity = np.array([0.0, 0.0, 0.0])
        self.positionOffset = np.array([0.0, 0.0, 0.0])
        
        self.angularVelocity = np.array([0.0, 0.0, 0.0])
        self.angularAcceleration = np.array([0.0, 0.0, 0.0])
        
        self.momentOfInertia = 1.0
        
        self.Force = None
        self.Torque = None
        
        self.xPoints=[]
        self.yPoints=[]
        self.zPoints=[]
        
        self.intensities = []
        self.forceVectors = []
        
        self.history = {
            "time": [],
            "position": [],
            "velocity": [],
            "force": [],
            "torque": [],
            "thetaX": [],
            "thetaY": [],
            "angularVelocity": []
        }
    def compute(self):
        self.xPoints=[]
        self.yPoints=[]
        self.zPoints=[]
        self.intensities = []
        self.forceVectors = []
        
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
                
                position += self.positionOffset
                
                r = np.sqrt(position[0]**2 + position[1]**2)
                I = gaussian(position[0], position[1])
                
                self.xPoints.append(position[0])
                self.yPoints.append(position[1])
                self.zPoints.append(position[2])
                
                self.intensities.append(I)
                
                normal = np.array([np.sin(self.thetaY) * np.cos(self.thetaX), -np.sin(self.thetaX), np.cos(self.thetaX) * np.cos(self.thetaY)])
                
                F = calcForce(I, normal, dA)
                
                self.forceVectors.append(F)
                
                tau = np.cross(position, F)
                
                Force = Force + F
                Torque = Torque + tau
        self.Force = Force
        self.Torque = Torque
        return Force, Torque

    def update(self, dt, time):
        acceleration = (self.Force / self.mass)
        self.angularAcceleration= (self.Torque / self.momentOfInertia)
        
        self.velocity += (acceleration * dt)
        self.angularVelocity += (self.angularAcceleration * dt)
        
        self.velocity *= 0.995 # Dampening, since tiny errors can accumulate and cause problems :3
        self.angularVelocity *= 0.995
        
        self.positionOffset += (self.velocity * dt)
        self.thetaX += (self.angularVelocity[0] * dt)
        self.thetaX += (self.angularVelocity[1] * dt)
        
        self.history["time"].append(time)
        self.history["position"].append(self.positionOffset.copy())
        self.history["velocity"].append(self.velocity.copy())
        self.history["force"].append(self.Force.copy())
        self.history["torque"].append(self.Torque.copy())
        self.history["thetaX"].append(self.thetaX)
        self.history["thetaY"].append(self.thetaY)
        self.history["angularVelocity"].append(self.angularVelocity.copy())

    def Visualize(self):
        fig = plt.figure()
        
        Force, Torque = self.compute()
        
        ax = fig.add_subplot(111, projection='3d')
        positions = np.array([self.history["position"]])
        
        ax.plot(positions[:,0], positions[:,1], positions[:,2])
        
        scatter = ax.scatter(self.xPoints, self.yPoints, self.zPoints, c=self.intensities, cmap='plasma', s=30)
        
        fig.colorbar(scatter, ax=ax, label='Beam Intensity (W/m^2)')
        
        ax.quiver(0, 0, 0, Force[0], Force[1], Force[2], color='red', length=2, normalize=True)
        
        ax.text2D(0.05, 0.95, f'Total Force: {Force}\nTotal Torque: {Torque}')
        
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
        
        self.mass = 1.0
        
        self.velocity = np.array([0.0, 0.0, 0.0])
        self.positionOffset = np.array([0.0, 0.0, 0.0])
        
        self.angularVelocity = np.array([0.0, 0.0, 0.0])
        self.angularAcceleration = np.array([0.0, 0.0, 0.0])
        
        self.momentOfInertia = 1.0
        
        self.Force = None
        self.Torque = None
        
        self.xPoints=[]
        self.yPoints=[]
        self.zPoints=[]
        
        self.intensities = []
        
        self.history = {
            "time": [],
            "position": [],
            "velocity": [],
            "force": [],
            "torque": [],
            "thetaX": [],
            "thetaY": [],
            "angularVelocity": []
        }

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
                
                position += self.positionOffset
                
                self.xPoints.append(position[0])
                self.yPoints.append(position[1])
                self.zPoints.append(position[2])
                
                normal = position / np.linalg.norm(position)
                
                r = np.sqrt(position[0]**2 + position[1]**2)
                I = gaussian(position[0], position[1])
                
                self.intensities.append(I)
                
                dtheta = (np.pi / 2) / self.resolution
                dphi = (2 * np.pi) / self.resolution
                
                dA = self.radius**2 * np.sin(theta) * dtheta * dphi
                
                F = calcForce(I, normal, dA)
                
                tau = np.cross(position, F)
                
                Force = Force + F
                Torque = Torque + tau
        self.Force = Force
        self.Torque = Torque
        return Force, Torque

    def update(self, dt, time):
        acceleration = (self.Force / self.mass)
        self.angularAcceleration= (self.Torque / self.momentOfInertia)
        
        self.velocity += (acceleration* dt)
        self.angularVelocity += (self.angularAcceleration * dt)
        
        self.velocity *= 0.995
        self.angularVelocity *= 0.995
        
        self.positionOffset += (self.velocity * dt)
        self.thetaX += (self.angularVelocity[0] * dt)
        self.thetaX += (self.angularVelocity[1] * dt)
        
        self.history["time"].append(time)
        self.history["position"].append(self.positionOffset.copy())
        self.history["velocity"].append(self.velocity.copy())
        self.history["force"].append(self.Force.copy())
        self.history["torque"].append(self.Torque.copy())
        self.history["thetaX"].append(self.thetaX)
        self.history["thetaY"].append(self.thetaY)
        self.history["angularVelocity"].append(self.angularVelocity.copy())
    
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
        
        self.mass = 1.0
        
        self.angularVelocity = np.array([0.0, 0.0, 0.0])
        self.angularAcceleration = np.array([0.0, 0.0, 0.0])
        
        self.momentOfInertia = 1.0
        
        self.velocity = np.array([0.0, 0.0, 0.0])
        self.positionOffset = np.array([0.0, 0.0, 0.0])
        
        self.Force = None
        self.Torque = None
        
        self.xPoints=[]
        self.yPoints=[]
        self.zPoints=[]
        self.intensities = []
        
        self.history = {
            "time": [],
            "position": [],
            "velocity": [],
            "force": [],    
            "torque": [],
            "thetaX": [],
            "thetaY": [],
            "angularVelocity": []
        }
    
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
                
                z = self.a * (self.radius**2 - r**2)
                
                normal = np.array([-2 * self.a * x, -2 * self.a * y, 1.0])
                normal = normal / np.linalg.norm(normal)
                
                normal = rotateX(normal, self.thetaX)
                normal = rotateY(normal, self.thetaY)
                
                position = np.array([x, y, z])
                
                position = rotateX(position, self.thetaX)
                position = rotateY(position, self.thetaY)
                
                position += self.positionOffset
                
                self.xPoints.append(position[0])
                self.yPoints.append(position[1])
                self.zPoints.append(position[2])
                
                I = gaussian(position[0], position[1])
                
                dr = self.radius / self.resolution
                dphi = 2 * np.pi / self.resolution
                
                dA = np.sqrt(1 + (2*self.a*position[0])**2 + (2*self.a*position[1])**2) * r * dr * dphi
                
                F = calcForce(I, normal, dA)
                
                tau = np.cross(position, F)
                
                Force += F
                Torque += tau
                
                self.intensities.append(I)
        self.Force=Force
        self.Torque=Torque
        return Force, Torque

    def update(self, dt, time):
        acceleration = (self.Force / self.mass)
        self.angularAcceleration= (self.Torque / self.momentOfInertia)
        
        self.velocity += (acceleration* dt)
        self.angularVelocity += (self.angularAcceleration * dt)
        
        self.velocity *= 0.995
        self.angularVelocity *= 0.995
        
        self.positionOffset += (self.velocity * dt)
        self.thetaX += (self.angularVelocity[0] * dt)
        self.thetaX += (self.angularVelocity[1] * dt)
        
        self.history["time"].append(time)
        self.history["position"].append(self.positionOffset.copy())
        self.history["velocity"].append(self.velocity.copy())
        self.history["force"].append(self.Force.copy())
        self.history["torque"].append(self.Torque.copy())
        self.history["thetaX"].append(self.thetaX)
        self.history["thetaY"].append(self.thetaY)
        self.history["angularVelocity"].append(self.angularVelocity.copy())

    def Visualize(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        scatter = ax.scatter(self.xPoints, self.yPoints, self.zPoints, c=self.intensities, cmap='plasma', s=30)
        
        fig.colorbar(scatter, ax=ax, label='Beam Intensity (W/m^2)')
        
        ax.quiver(0, 0, 0, self.Force[0], self.Force[1], self.Force[2], color='red',length=2, normalize=True)
        
        ax.set_box_aspect([1,1,1])
        ax.view_init(elev=30, azim=45)
        
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        
        ax.margins(0.5)
        
        plt.show()

def simLoop(timestep, steps, sail):
    time = 0
    dt = timestep
    print("="*90)
    for step in range(steps):
        # print(f"{sail.positionOffset} (Time: {time})")
        
        sail.compute()
        sail.update(dt, time)
        
        time += dt
    sail.Visualize()



nyxRect1 = RectangleLightSail(10, 10, 100, 0, 0, np.radians(0), np.radians(0))
nyxRect2 = RectangleLightSail(10, 10, 100, 0, 0, np.radians(20), np.radians(0))
nyxRect3 = RectangleLightSail(10, 10, 100, 0, 0, np.radians(45), np.radians(0))
nyxRect4 = RectangleLightSail(10, 10, 100, 0, 0, np.radians(90), np.radians(0))

nyxSphere1 = SphereLightSail(1, 100, np.radians(0), np.radians(0))
nyxSphere2 = SphereLightSail(1, 100, np.radians(30), np.radians(0))
nyxSphere3 = SphereLightSail(1, 100, np.radians(60), np.radians(0))

nyxParaboloid1 = ParaboloidLightSail(1, 0.5, 100, np.radians(0), np.radians(0))

nyxDisc1 = ParaboloidLightSail(5, 0, 100, np.radians(0), np.radians(0))

dt = 0.01
simLoop(dt, 50, nyxRect1)
# simLoop(dt, 50, nyxSphere1)
# simLoop(dt, 50, nyxParaboloid1)
# simLoop(dt, 50, nyxDisc1)

# TotalForce, TotalTorque = nyxParaboloid1.compute()
# print(f"Paraboloid Reflector Test Case 1: \n ================== \n Total force is {TotalForce} \n Total torque is {TotalTorque} \n ==================")
# nyxParaboloid1.Visualize()

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