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
        F = 2 * (I / c) * (cosTheta ** 2) * dA * normal
    else:
        F = np.array([0.0, 0.0, 0.0])
    
    return F


def _is_centered_untilted(sail, tol=1e-8):
    centered = (
        abs(getattr(sail, 'thetaX', 0.0)) < tol
        and abs(getattr(sail, 'thetaY', 0.0)) < tol
        and abs(getattr(sail, 'positionOffset', np.array([0.0,0.0,0.0]))[0]) < tol
        and abs(getattr(sail, 'positionOffset', np.array([0.0,0.0,0.0]))[1]) < tol
    )
    if hasattr(sail, 'CenterX') and hasattr(sail, 'CenterY'):
        centered = (
            centered
            and abs(getattr(sail, 'CenterX', 0.0)) < tol
            and abs(getattr(sail, 'CenterY', 0.0)) < tol
        )
    return centered

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
        
        self.momentOfInertia = (1/12) * self.mass * (self.width**2 + self.height**2)
        
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
                
                rotatedPos = position.copy()
                
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
                
                tau = np.cross(rotatedPos, F) 
                
                Force = Force + F
                Torque = Torque + tau
        self.Force = Force
        self.Torque = Torque

        if _is_centered_untilted(self):
            self.Force[0] = 0.0
            self.Force[1] = 0.0
            self.Torque[0] = 0.0
            self.Torque[1] = 0.0
        else:
            tol = 1e-14
            self.Force[np.abs(self.Force) < tol] = 0.0
            self.Torque[np.abs(self.Torque) < tol] = 0.0

        return self.Force, self.Torque

    def derivatives(self, position, velocity, thetaX, thetaY, angularVelocity):
        oldPosition = self.positionOffset.copy()
        oldVelocity = self.velocity.copy()
        oldThetaX = self.thetaX
        oldThetaY = self.thetaY
        oldAngVel = self.angularVelocity.copy()
        
        self.positionOffset = position
        self.velocity = velocity
        self.thetaX = thetaX
        self.thetaY = thetaY
        self.angularVelocity = angularVelocity
        
        self.compute()
        
        dpos = self.velocity.copy() #dx/dt = v
        dvel = self.Force / self.mass #dv/dt = a
        dthX = self.angularVelocity[0] #d(thetaX)/dt = omegaX
        dthY = self.angularVelocity[1] #d(thetaY)/dt = omegaY
        dangVel = self.Torque / self.momentOfInertia #d(omega)/dt = alpha # lmao "dang"
        
        self.positionOffset = oldPosition
        self.velocity = oldVelocity
        self.thetaX = oldThetaX
        self.thetaY = oldThetaY
        self.angularVelocity = oldAngVel
        
        return dpos, dvel, dthX, dthY, dangVel
    
    def update(self, dt, time):
        pos = self.positionOffset.copy()
        velocity = self.velocity.copy()
        thetaX = self.thetaX
        thetaY = self.thetaY
        angularVelocity = self.angularVelocity.copy()
        
        #k1
        dpos1, dvel1, dthX1, dthY1, dangVel1 = self.derivatives(pos, velocity, thetaX, thetaY, angularVelocity)
        
        #k2
        p2 = pos + 0.5 * dt * dpos1
        v2 = velocity + 0.5 * dt * dvel1
        thx2 = thetaX + 0.5 * dt * dthX1
        thy2 = thetaY + 0.5 * dt * dthY1
        w2 = angularVelocity + 0.5 * dt * dangVel1
        
        dpos2, dvel2, dthX2, dthY2, dangVel2 = self.derivatives(p2, v2, thx2, thy2, w2)
        
        #k3
        p3 = pos + 0.5 * dt * dpos2
        v3 = velocity + 0.5 * dt * dvel2
        thx3 = thetaX + 0.5 * dt * dthX2
        thy3 = thetaY + 0.5 * dt * dthY2
        w3 = angularVelocity + 0.5 * dt * dangVel2
        
        dpos3, dvel3, dthX3, dthY3, dangVel3 = self.derivatives(p3, v3, thx3, thy3, w3)
        
        #k4
        p4 = pos + dt * dpos3
        v4 = velocity + dt * dvel3
        thx4 = thetaX + dt * dthX3
        thy4 = thetaY + dt * dthY3
        w4 = angularVelocity + dt * dangVel3
        
        dpos4, dvel4, dthX4, dthY4, dangVel4 = self.derivatives(p4, v4, thx4, thy4, w4)
        
        #weighted average or something idk im just following a random RK4 tutorial i aint that smart
        finalpos = pos + (dt/ 6.0) * (dpos1 + 2*dpos2 + 2*dpos3 + dpos4)
        finalvelocity = velocity + (dt/ 6.0) * (dvel1 + 2*dvel2 + 2*dvel3 + dvel4)
        finalthetaX = thetaX + (dt/ 6.0) * (dthX1 + 2*dthX2 + 2*dthX3 + dthX4)
        finalthetaY = thetaY + (dt/ 6.0) * (dthY1 + 2*dthY2 + 2*dthY3 + dthY4)
        finalangVel = angularVelocity + (dt/ 6.0) * (dangVel1 + 2*dangVel2 + 2*dangVel3 + dangVel4)
        
        self.positionOffset = finalpos
        self.velocity = finalvelocity
        self.thetaX = finalthetaX
        self.thetaY = finalthetaY
        self.angularVelocity = finalangVel

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
        positions = np.array(self.history["position"])
        
        ax.plot(positions[:,0], positions[:,1], positions[:,2])
        
        scatter = ax.scatter(self.xPoints, self.yPoints, self.zPoints, c=self.intensities, cmap='plasma', s=30)
        
        fig.colorbar(scatter, ax=ax, label='Beam Intensity (W/m^2)')
        
        ax.quiver(0, 0, 0, self.Force[0], self.Force[1], self.Force[2], color='red', length=2, normalize=True)
        
        ax.text2D(0.05, 0.95, f'Total Force: {self.Force}\nTotal Torque: {self.Torque}')
        
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
        self.intensities = []
        
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
        
        if _is_centered_untilted(self):
            self.Force[0] = 0.0
            self.Force[1] = 0.0
            self.Torque[0] = 0.0
            self.Torque[1] = 0.0
        else:
            tol = 1e-14
            self.Force[np.abs(self.Force) < tol] = 0.0
            self.Torque[np.abs(self.Torque) < tol] = 0.0
        return self.Force, self.Torque

    def derivatives(self, position, velocity, thetaX, thetaY, angularVelocity):
        oldPosition = self.positionOffset.copy()
        oldVelocity = self.velocity.copy()
        oldThetaX = self.thetaX
        oldThetaY = self.thetaY
        oldAngVel = self.angularVelocity.copy()
        
        self.positionOffset = position
        self.velocity = velocity
        self.thetaX = thetaX
        self.thetaY = thetaY
        self.angularVelocity = angularVelocity
        
        self.compute()
        
        dpos = self.velocity.copy() #dx/dt = v
        dvel = self.Force / self.mass #dv/dt = a
        dthX = self.angularVelocity[0] #d(thetaX)/dt = omegaX
        dthY = self.angularVelocity[1] #d(thetaY)/dt = omegaY
        dangVel = self.Torque / self.momentOfInertia #d(omega)/dt = alpha # lmao "dang"
        
        self.positionOffset = oldPosition
        self.velocity = oldVelocity
        self.thetaX = oldThetaX
        self.thetaY = oldThetaY
        self.angularVelocity = oldAngVel
        
        return dpos, dvel, dthX, dthY, dangVel
    
    def update(self, dt, time):
        pos = self.positionOffset.copy()
        velocity = self.velocity.copy()
        thetaX = self.thetaX
        thetaY = self.thetaY
        angularVelocity = self.angularVelocity.copy()
        
        #k1
        dpos1, dvel1, dthX1, dthY1, dangVel1 = self.derivatives(pos, velocity, thetaX, thetaY, angularVelocity)
        
        #k2
        p2 = pos + 0.5 * dt * dpos1
        v2 = velocity + 0.5 * dt * dvel1
        thx2 = thetaX + 0.5 * dt * dthX1
        thy2 = thetaY + 0.5 * dt * dthY1
        w2 = angularVelocity + 0.5 * dt * dangVel1
        
        dpos2, dvel2, dthX2, dthY2, dangVel2 = self.derivatives(p2, v2, thx2, thy2, w2)
        
        #k3
        p3 = pos + 0.5 * dt * dpos2
        v3 = velocity + 0.5 * dt * dvel2
        thx3 = thetaX + 0.5 * dt * dthX2
        thy3 = thetaY + 0.5 * dt * dthY2
        w3 = angularVelocity + 0.5 * dt * dangVel2
        
        dpos3, dvel3, dthX3, dthY3, dangVel3 = self.derivatives(p3, v3, thx3, thy3, w3)
        
        #k4
        p4 = pos + dt * dpos3
        v4 = velocity + dt * dvel3
        thx4 = thetaX + dt * dthX3
        thy4 = thetaY + dt * dthY3
        w4 = angularVelocity + dt * dangVel3
        
        dpos4, dvel4, dthX4, dthY4, dangVel4 = self.derivatives(p4, v4, thx4, thy4, w4)
        
        #weighted average or something idk im just following a random RK4 tutorial i aint that smart
        finalpos = pos + (dt/ 6.0) * (dpos1 + 2*dpos2 + 2*dpos3 + dpos4)
        finalvelocity = velocity + (dt/ 6.0) * (dvel1 + 2*dvel2 + 2*dvel3 + dvel4)
        finalthetaX = thetaX + (dt/ 6.0) * (dthX1 + 2*dthX2 + 2*dthX3 + dthX4)
        finalthetaY = thetaY + (dt/ 6.0) * (dthY1 + 2*dthY2 + 2*dthY3 + dthY4)
        finalangVel = angularVelocity + (dt/ 6.0) * (dangVel1 + 2*dangVel2 + 2*dangVel3 + dangVel4)
        
        self.positionOffset = finalpos
        self.velocity = finalvelocity
        self.thetaX = finalthetaX
        self.thetaY = finalthetaY
        self.angularVelocity = finalangVel

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
        
        positions = np.array(self.history["position"])
        
        ax.plot(positions[:,0], positions[:,1], positions[:,2])
        
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
                
                z = -self.a * (self.radius**2 - r**2)
                
                normal = np.array([2 * self.a * x, 2 * self.a * y, 1.0])
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
                
                dA = np.sqrt(1 + (2*self.a*x)**2 + (2*self.a*y)**2) * r * dr * dphi
                
                F = calcForce(I, normal, dA)
                
                tau = np.cross(position, F)
                
                Force += F
                Torque += tau
                
                self.intensities.append(I)
        self.Force = Force
        self.Torque = Torque
        
        if _is_centered_untilted(self):
            self.Force[0] = 0.0
            self.Force[1] = 0.0
            self.Torque[0] = 0.0
            self.Torque[1] = 0.0
        else:
            tol = 1e-14
            self.Force[np.abs(self.Force) < tol] = 0.0
            self.Torque[np.abs(self.Torque) < tol] = 0.0
        return self.Force, self.Torque

    def derivatives(self, position, velocity, thetaX, thetaY, angularVelocity):
        oldPosition = self.positionOffset.copy()
        oldVelocity = self.velocity.copy()
        oldThetaX = self.thetaX
        oldThetaY = self.thetaY
        oldAngVel = self.angularVelocity.copy()
        
        self.positionOffset = position
        self.velocity = velocity
        self.thetaX = thetaX
        self.thetaY = thetaY
        self.angularVelocity = angularVelocity
        
        self.compute()
        
        dpos = self.velocity.copy() #dx/dt = v
        dvel = self.Force / self.mass #dv/dt = a
        dthX = self.angularVelocity[0] #d(thetaX)/dt = omegaX
        dthY = self.angularVelocity[1] #d(thetaY)/dt = omegaY
        dangVel = self.Torque / self.momentOfInertia #d(omega)/dt = alpha # lmao "dang"
        
        self.positionOffset = oldPosition
        self.velocity = oldVelocity
        self.thetaX = oldThetaX
        self.thetaY = oldThetaY
        self.angularVelocity = oldAngVel
        
        return dpos, dvel, dthX, dthY, dangVel
    
    def update(self, dt, time):
        pos = self.positionOffset.copy()
        velocity = self.velocity.copy()
        thetaX = self.thetaX
        thetaY = self.thetaY
        angularVelocity = self.angularVelocity.copy()
        
        #k1
        dpos1, dvel1, dthX1, dthY1, dangVel1 = self.derivatives(pos, velocity, thetaX, thetaY, angularVelocity)
        
        #k2
        p2 = pos + 0.5 * dt * dpos1
        v2 = velocity + 0.5 * dt * dvel1
        thx2 = thetaX + 0.5 * dt * dthX1
        thy2 = thetaY + 0.5 * dt * dthY1
        w2 = angularVelocity + 0.5 * dt * dangVel1
        
        dpos2, dvel2, dthX2, dthY2, dangVel2 = self.derivatives(p2, v2, thx2, thy2, w2)
        
        #k3
        p3 = pos + 0.5 * dt * dpos2
        v3 = velocity + 0.5 * dt * dvel2
        thx3 = thetaX + 0.5 * dt * dthX2
        thy3 = thetaY + 0.5 * dt * dthY2
        w3 = angularVelocity + 0.5 * dt * dangVel2
        
        dpos3, dvel3, dthX3, dthY3, dangVel3 = self.derivatives(p3, v3, thx3, thy3, w3)
        
        #k4
        p4 = pos + dt * dpos3
        v4 = velocity + dt * dvel3
        thx4 = thetaX + dt * dthX3
        thy4 = thetaY + dt * dthY3
        w4 = angularVelocity + dt * dangVel3
        
        dpos4, dvel4, dthX4, dthY4, dangVel4 = self.derivatives(p4, v4, thx4, thy4, w4)
        
        #weighted average or something idk im just following a random RK4 tutorial i aint that smart
        finalpos = pos + (dt/ 6.0) * (dpos1 + 2*dpos2 + 2*dpos3 + dpos4)
        finalvelocity = velocity + (dt/ 6.0) * (dvel1 + 2*dvel2 + 2*dvel3 + dvel4)
        finalthetaX = thetaX + (dt/ 6.0) * (dthX1 + 2*dthX2 + 2*dthX3 + dthX4)
        finalthetaY = thetaY + (dt/ 6.0) * (dthY1 + 2*dthY2 + 2*dthY3 + dthY4)
        finalangVel = angularVelocity + (dt/ 6.0) * (dangVel1 + 2*dangVel2 + 2*dangVel3 + dangVel4)
        
        self.positionOffset = finalpos
        self.velocity = finalvelocity
        self.thetaX = finalthetaX
        self.thetaY = finalthetaY
        self.angularVelocity = finalangVel

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
        
        positions = np.array(self.history["position"])
        
        ax.plot(positions[:,0], positions[:,1], positions[:,2])
        
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
        print(f"{sail.positionOffset} (Time: {time})")
        sail.compute()
        sail.update(dt, time)
        
        time += dt
    sail.Visualize()


if __name__ == '__main__':

    nyxRect1 = RectangleLightSail(10, 10, 100, 0, 0, np.radians(0), np.radians(0))
    nyxRect2 = RectangleLightSail(10, 10, 100, 0, 0, np.radians(20), np.radians(0))
    nyxRect3 = RectangleLightSail(10, 10, 100, 0, 0, np.radians(45), np.radians(0))
    nyxRect4 = RectangleLightSail(10, 10, 100, 0, 0, np.radians(90), np.radians(0))

    nyxSphere1 = SphereLightSail(1, 100, np.radians(0), np.radians(0))
    nyxSphere2 = SphereLightSail(1, 100, np.radians(30), np.radians(0))
    nyxSphere3 = SphereLightSail(1, 100, np.radians(60), np.radians(0))

    nyxParaboloid1 = ParaboloidLightSail(1, 0.5, 100, np.radians(0), np.radians(0))

    nyxDisc1 = ParaboloidLightSail(5, 0, 100, np.radians(0), np.radians(0))

    dt = 1.0
    steps = 100
    simLoop(dt, steps, nyxRect1)
    simLoop(dt, steps, nyxParaboloid1)
    # simLoop(dt, steps, nyxSphere1)
    # simLoop(dt, steps, nyxDisc1)

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