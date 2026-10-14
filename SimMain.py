import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal

from materials import Graphene, Aluminium, SiliconCarbide
from analysis.plotting import *

# NOTE FOR SELF: Run using & "C:\Users\dell\AppData\Local\Programs\Python\Python313\python.exe" "C:\Users\dell\Documents\Makeiums\Nyx\Simulation\SimMain.py"

#Constants
c = 3e8 #The Speed of Light in m/s. This is a constant for the simulation, and can be changed to test different scenarios.

I0 = 1e9 #The Intensity at The Exact Center of the Beam. This is a constant for the simulation, and can be changed to test different scenarios.
w0 = 1.0 #The beam waist, which is the radius at which the intensity falls to 1/e^2 of its maximum value. This is also a constant for the simulation, and can be changed to test different scenarios.
lambdaLaser = 1e-6 # 1 micron laser
zR = np.pi * w0**2/lambdaLaser #Rayleigh range
sigma = 5.670374419e-8

beamX = 0
beamY = 0
beamDih = np.array([0.0, 0.0, 1.0]) #Beam direction, i thought it would be funny to name it "Beam Dih" 
beamDih = beamDih / np.linalg.norm(beamDih)


def beamradi(z):
    return w0 * np.sqrt(1 + (z / zR)**2)

#Gaussian Beam Equation: I(r) = I0 * exp(-2 * r^2 / w0^2), How the intensity of the beam changes with distance from the center. This is a fundamental equation for simulating the interaction of light with the sail, and is used to calculate the force and torque on each element of the sail based on its position relative to the center of the beam.
def gaussian(x, y, z=0):
    wz = beamradi(z)
    
    r = np.sqrt((x - beamX)**2 + (y - beamY)**2)
    
    return I0 * (w0 / wz)**2 * np.exp(-2 * r**2 / wz**2)

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

def quaternionToRotationMatrix(q):
    q = normalizeQuaternion(q)
    
    w, x, y, z = q
    
    return np.array([
        [1 - 2*(y**2 + z**2), 2*(x*y - z*w), 2*(x*z + y*w)],
        [2*(x*y + z*w), 1 - 2*(x**2 + z**2), 2*(y*z - x*w)],
        [2*(x*z - y*w), 2*(y*z + x*w), 1 - 2*(x**2 + y**2)]
    ])

def calcForce(I, normal, dA, reflectivity, absorptivity):
    if not np.all(np.isfinite(normal)):
        return np.array([0.0, 0.0, 0.0])
    
    norm = np.linalg.norm(normal)
    
    if not np.isfinite(norm) or norm < 1e-15:
        return np.array([0.0, 0.0, 0.0])
    normal = normal / norm

    cosTheta = np.dot(normal, beamDih)
    
    incoming = -beamDih / np.linalg.norm(beamDih)
    
    reflected = incoming - 2 * np.dot(incoming, normal) * normal
    reflected =  reflected / np.linalg.norm(reflected)
    
    if not np.isfinite(cosTheta):
        return np.array([0.0, 0.0, 0.0])
    
    
    cosTheta = np.clip(cosTheta, -1.0, 1.0)

    if cosTheta <= 1e-10:
        return np.array([0.0, 0.0, 0.0])
    
    deltaP = reflected - incoming
    fmag = (I / c) * cosTheta * dA
    
    reflectedForce = reflectivity * fmag * deltaP
    absorbedForce = absorptivity * fmag * beamDih
    
    return reflectedForce + absorbedForce


def _is_centered_untilted(sail, tol=1e-8):
    pos = getattr(sail, 'positionOffset', np.array([0.0, 0.0, 0.0]))
    centered = (
        abs(pos[0]) < tol
        and abs(pos[1]) < tol
        and abs(beamX) < tol
        and abs(beamY) < tol
    )
    if hasattr(sail, 'CenterX') and hasattr(sail, 'CenterY'):
        centered = (
            centered
            and abs(getattr(sail, 'CenterX', 0.0)) < tol
            and abs(getattr(sail, 'CenterY', 0.0)) < tol
        )

    orientation = getattr(sail, 'orientation', np.array([1.0, 0.0, 0.0, 0.0]))
    orientation = normalizeQuaternion(orientation)
    untilted = np.allclose(orientation, np.array([1.0, 0.0, 0.0, 0.0]), atol=tol)
    return centered and untilted

def quaternionMultiply(q1, q2):
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    
    return np.array([
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2
    ])

def normalizeQuaternion(q):
    norm = np.linalg.norm(q)
    if not np.isfinite(norm) or norm < 1e-15:
        return np.array([1.0, 0.0, 0.0, 0.0])
    return q / norm

def rotateVector(vector, quaternion):
    q = quaternion
    if not np.all(np.isfinite(q)):
        q = np.array([1.0, 0.0, 0.0, 0.0])
    q = normalizeQuaternion(q)
    
    qConj = np.array([q[0], -q[1], -q[2], -q[3]])
    
    vecQuat = np.array([0.0, vector[0], vector[1], vector[2]])
    
    rotated = quaternionMultiply(quaternionMultiply(q, vecQuat), qConj)
    
    return rotated[1:]

def angularVelocityToQuaternion(omega, dt):
    angle = np.linalg.norm(omega) * dt
    
    if angle < 1e-8:
        return np.array([1.0, 0.0, 0.0, 0.0])
    
    axis = omega / np.linalg.norm(omega)
    sin_half = np.sin(angle / 2)
    return np.array([np.cos(angle / 2), axis[0] * sin_half, axis[1] * sin_half, axis[2] * sin_half])

def eul2quat(thetaX, thetaY, thetaZ=0):
    cx = np.cos(thetaX/2)
    sx = np.sin(thetaX/2)
    
    cy = np.cos(thetaY/2)
    sy = np.sin(thetaY/2)
    
    cz = np.cos(thetaZ/2)
    sz = np.sin(thetaZ/2)
    
    w = cx*cy*cz + sx*sy*sz
    x = sx*cy*cz - cx*sy*sz
    y = cx*sy*cz + sx*cy*sz
    z = cx*cy*sz - sx*sy*cz
    
    return np.array([w, x, y, z])


class RectangleLightSail():
    def __init__(self, width, height, resolution, StartingX=0, StartingY=0, thetaX=0, thetaY=0, thickness=1e-6, material=Graphene):
        self.width = width
        self.height = height
        self.resolution = resolution
        self.CenterX = StartingX
        self.CenterY = StartingY
        self.thickness = thickness
        self.material = material
        self.mass = self.material.density * self.width * self.height * self.thickness
        
        self.velocity = np.array([0.0, 0.0, 0.0])
        self.positionOffset = np.array([0.0, 0.0, 0.0])
        
        self.angularVelocity = np.array([0.0, 0.0, 0.0])
        self.angularAcceleration = np.array([0.0, 0.0, 0.0])
        
        Ixx = (1/12) * self.mass * (self.height**2)
        Iyy = (1/12) * self.mass * (self.width**2)
        Izz = (1/12) * self.mass * (self.width**2 + self.height**2)
        
        self.inertiaTensor = np.array([
            [Ixx, 0.0, 0.0],
            [0.0, Iyy, 0.0],
            [0.0, 0.0, Izz]
        ])
        
        
        self.reflectivity = self.material.reflectivity
        self.absorptivity = self.material.absorptivity
        
        self.thetaX = thetaX
        self.thetaY = thetaY
        self.thetaZ = 0.0
        self.orientation = eul2quat(thetaX, thetaY) #w, x, y, z
        
        self.Force = None
        self.Torque = None
        
        self.xPoints=[]
        self.yPoints=[]
        self.zPoints=[]
        
        self.intensities = []
        self.forceVectors = []
        
        self.temperatureMap = []
        
        self.history = {
            "time": [],
            "position": [],
            "velocity": [],
            "force": [],
            "torque": [],
            "orientation": [],
            "angularVelocity": [],
            "maxTemperature": [],
            "averageTemperature": []
        }
    def compute(self):
        self.xPoints=[]
        self.yPoints=[]
        self.zPoints=[]
        self.intensities = []
        self.forceVectors = []
        self.temperatureMap = []
        
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
                normal = np.array([0.0, 0.0, 1.0])
                
                position = rotateVector(position, self.orientation)
                
                rotatedPos = position.copy()
                
                position += self.positionOffset
                
                r = np.sqrt(position[0]**2 + position[1]**2)
                I = gaussian(position[0], position[1], position[2])
                absorbedFlux = self.absorptivity * I
                temp = (absorbedFlux / (self.material.emissivity * sigma)) ** 0.25
                
                self.xPoints.append(position[0])
                self.yPoints.append(position[1])
                self.zPoints.append(position[2])
                
                self.intensities.append(I)
                self.temperatureMap.append(temp)
                
                normal = rotateVector(normal, self.orientation)
                
                F = calcForce(I, normal, dA, self.reflectivity, self.absorptivity)
                
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

    def derivatives(self, position, velocity, orientation, angularVelocity):
        oldPosition = self.positionOffset.copy()
        oldVelocity = self.velocity.copy()
        oldOrientation = self.orientation.copy()
        oldAngVel = self.angularVelocity.copy()
        
        self.positionOffset = position
        self.velocity = velocity
        self.orientation = orientation
        self.angularVelocity = angularVelocity
        
        self.compute()
        
        dpos = self.velocity.copy() #dx/dt = v
        dvel = self.Force / self.mass #dv/dt = a
        omegaQuat = np.array([0.0, angularVelocity[0], angularVelocity[1], angularVelocity[2]])
        qDot = 0.5 * quaternionMultiply(omegaQuat, orientation)
        # dangVel = np.linalg.solve(self.inertiaTensor, self.Torque) #d(omega)/dt = alpha # lmao "dang"
        
        R = quaternionToRotationMatrix(orientation)
        
        I_world = R @ self.inertiaTensor @ R.T
        angularMomentum = I_world @ angularVelocity
        
        dangVel = np.linalg.solve(I_world, self.Torque - np.cross(angularVelocity, angularMomentum))
        
        self.positionOffset = oldPosition
        self.velocity = oldVelocity
        self.orientation = oldOrientation
        self.angularVelocity = oldAngVel
        
        return dpos, dvel, qDot, dangVel
    
    def update(self, dt, time):
        pos = self.positionOffset.copy()
        velocity = self.velocity.copy()
        orientation = self.orientation.copy()
        angularVelocity = self.angularVelocity.copy()
        
        #k1
        dpos1, dvel1, dq1, dangVel1 = self.derivatives(pos, velocity, orientation, angularVelocity)
        
        #k2
        p2 = pos + 0.5 * dt * dpos1
        v2 = velocity + 0.5 * dt * dvel1
        q2 = orientation + 0.5 * dt * dq1
        q2 = normalizeQuaternion(q2)
        w2 = angularVelocity + 0.5 * dt * dangVel1
        
        dpos2, dvel2, dq2, dangVel2 = self.derivatives(p2, v2, q2, w2)
        
        #k3
        p3 = pos + 0.5 * dt * dpos2
        v3 = velocity + 0.5 * dt * dvel2
        q3 = orientation + 0.5 * dt * dq2
        q3 = normalizeQuaternion(q3)
        w3 = angularVelocity + 0.5 * dt * dangVel2
        
        dpos3, dvel3, dq3, dangVel3 = self.derivatives(p3, v3, q3, w3)
        
        #k4
        p4 = pos + dt * dpos3
        v4 = velocity + dt * dvel3
        q4 = orientation + dt * dq3
        q4 = normalizeQuaternion(q4)
        
        w4 = angularVelocity + dt * dangVel3
        
        dpos4, dvel4, dq4, dangVel4 = self.derivatives(p4, v4, q4, w4)
        
        #weighted average or something idk im just following a random RK4 tutorial i aint that smart
        finalpos = pos + (dt/ 6.0) * (dpos1 + 2*dpos2 + 2*dpos3 + dpos4)
        finalvelocity = velocity + (dt/ 6.0) * (dvel1 + 2*dvel2 + 2*dvel3 + dvel4)
        finalOrientation = orientation + (dt / 6.0) * (dq1 + 2*dq2 + 2*dq3 + dq4)
        finalOrientation = normalizeQuaternion(finalOrientation)
        finalangVel = angularVelocity + (dt/ 6.0) * (dangVel1 + 2*dangVel2 + 2*dangVel3 + dangVel4)
        
        self.positionOffset = finalpos
        self.velocity = finalvelocity
        self.orientation = finalOrientation
        self.angularVelocity = finalangVel

        if self.temperatureMap:
            self.history["maxTemperature"].append(float(np.max(self.temperatureMap)))
            self.history["averageTemperature"].append(float(np.mean(self.temperatureMap)))
        else:
            self.history["maxTemperature"].append(0.0)
            self.history["averageTemperature"].append(0.0)
        
        self.history["time"].append(float(Decimal(str(time)) + Decimal(str(dt))))
        self.history["position"].append(self.positionOffset.copy())
        self.history["velocity"].append(self.velocity.copy())
        self.history["force"].append(self.Force.copy())
        self.history["torque"].append(self.Torque.copy())
        self.history["orientation"].append(self.orientation.copy())

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
    def __init__(self, radius, resolution, thetaX, thetaY, thickness=1e-6, material = Graphene):
        self.radius = radius
        self.resolution = resolution
        self.thickness = thickness
        self.material = material
        
        self.area = 2 * np.pi * self.radius**2
        self.mass = self.material.density * self.area * self.thickness
        
        self.velocity = np.array([0.0, 0.0, 0.0])
        self.positionOffset = np.array([0.0, 0.0, 0.0])
        
        self.angularVelocity = np.array([0.0, 0.0, 0.0])
        self.angularAcceleration = np.array([0.0, 0.0, 0.0])
        
        Ixx = (5/12) * self.mass * (self.radius**2)
        Izz = (2/3) * self.mass * (self.radius**2)
        
        self.inertiaTensor = np.array([
            [Ixx, 0.0, 0.0],
            [0.0, Ixx, 0.0],
            [0.0, 0.0, Izz]
        ])
        
        self.reflectivity = self.material.reflectivity
        self.absorptivity = self.material.absorptivity
        self.orientation = eul2quat(thetaX, thetaY) #w, x, y, z
        
        self.Force = None
        self.Torque = None
        
        self.xPoints=[]
        self.yPoints=[]
        self.zPoints=[]
        
        self.intensities = []
        self.temperatureMap = []
        
        self.history = {
            "time": [],
            "position": [],
            "velocity": [],
            "force": [],
            "torque": [],
            "orientation": [],
            "angularVelocity": [],
            "maxTemperature": [],
            "averageTemperature": []
        }

    def compute(self):
        self.xPoints=[]
        self.yPoints=[]
        self.zPoints=[]
        self.intensities = []
        self.temperatureMap = []
        
        Force = np.array([0.0, 0.0, 0.0])
        Torque = np.array([0.0, 0.0, 0.0])
        
        for i in range(self.resolution):
            for j in range(self.resolution):
                
                theta = (np.pi / 2) * (i + 0.5) / self.resolution
                phi = 2 * np.pi * (j + 0.5) / self.resolution
                
                x = self.radius * np.sin(theta) * np.cos(phi)
                y = self.radius * np.sin(theta) * np.sin(phi)
                z = self.radius * np.cos(theta)
                
                localPosition = np.array([x, y, z - self.radius / 2])
                
                normal = np.array([x, y, z]) / self.radius
                
                rotatedPos = rotateVector(localPosition, self.orientation)
                normal = rotateVector(normal, self.orientation)
                
                position = rotatedPos + self.positionOffset
                
                self.xPoints.append(position[0])
                self.yPoints.append(position[1])
                self.zPoints.append(position[2])
                
                r = np.sqrt(position[0]**2 + position[1]**2)
                I = gaussian(position[0], position[1], position[2])
                self.intensities.append(I)
                
                dtheta = (np.pi / 2) / self.resolution
                dphi = (2 * np.pi) / self.resolution
                
                dA = self.radius**2 * np.sin(theta) * dtheta * dphi
                
                absorbedFlux = self.absorptivity * I
                temp = (absorbedFlux / (self.material.emissivity * sigma)) ** 0.25
                
                self.temperatureMap.append(temp)
                
                F = calcForce(I, normal, dA, self.reflectivity, self.absorptivity)
                
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

    def derivatives(self, position, velocity, orientation, angularVelocity):
        oldPosition = self.positionOffset.copy()
        oldVelocity = self.velocity.copy()
        oldOrientation = self.orientation.copy()
        oldAngVel = self.angularVelocity.copy()
        
        self.positionOffset = position
        self.velocity = velocity
        self.orientation = orientation
        self.angularVelocity = angularVelocity
        
        self.compute()
        
        dpos = self.velocity.copy() #dx/dt = v
        dvel = self.Force / self.mass #dv/dt = a
        omegaQuat = np.array([0.0, angularVelocity[0], angularVelocity[1], angularVelocity[2]])
        qDot = 0.5 * quaternionMultiply(omegaQuat, orientation)
        # dangVel = np.linalg.solve(self.inertiaTensor, self.Torque) #d(omega)/dt = alpha # lmao "dang"
        
        R = quaternionToRotationMatrix(orientation)
        
        I_world = R @ self.inertiaTensor @ R.T
        angularMomentum = I_world @ angularVelocity
        
        dangVel = np.linalg.solve(I_world, self.Torque - np.cross(angularVelocity, angularMomentum))
        
        self.positionOffset = oldPosition
        self.velocity = oldVelocity
        self.orientation = oldOrientation
        self.angularVelocity = oldAngVel
        
        return dpos, dvel, qDot, dangVel
    
    def update(self, dt, time):
        pos = self.positionOffset.copy()
        velocity = self.velocity.copy()
        orientation = self.orientation.copy()
        angularVelocity = self.angularVelocity.copy()
        
        #k1
        dpos1, dvel1, dq1, dangVel1 = self.derivatives(pos, velocity, orientation, angularVelocity)
        
        #k2
        p2 = pos + 0.5 * dt * dpos1
        v2 = velocity + 0.5 * dt * dvel1
        q2 = orientation + 0.5 * dt * dq1
        q2 = normalizeQuaternion(q2)
        w2 = angularVelocity + 0.5 * dt * dangVel1
        
        dpos2, dvel2, dq2, dangVel2 = self.derivatives(p2, v2, q2, w2)
        
        #k3
        p3 = pos + 0.5 * dt * dpos2
        v3 = velocity + 0.5 * dt * dvel2
        q3 = orientation + 0.5 * dt * dq2
        q3 = normalizeQuaternion(q3)
        w3 = angularVelocity + 0.5 * dt * dangVel2
        
        dpos3, dvel3, dq3, dangVel3 = self.derivatives(p3, v3, q3, w3)
        
        #k4
        p4 = pos + dt * dpos3
        v4 = velocity + dt * dvel3
        q4 = orientation + dt * dq3
        q4 = normalizeQuaternion(q4)
        q4 = normalizeQuaternion(q4)
        
        w4 = angularVelocity + dt * dangVel3
        
        dpos4, dvel4, dq4, dangVel4 = self.derivatives(p4, v4, q4, w4)
        
        #weighted average or something idk im just following a random RK4 tutorial i aint that smart
        finalpos = pos + (dt/ 6.0) * (dpos1 + 2*dpos2 + 2*dpos3 + dpos4)
        finalvelocity = velocity + (dt/ 6.0) * (dvel1 + 2*dvel2 + 2*dvel3 + dvel4)
        finalOrientation = orientation + (dt / 6.0) * (dq1 + 2*dq2 + 2*dq3 + dq4)
        finalOrientation = normalizeQuaternion(finalOrientation)
        finalangVel = angularVelocity + (dt/ 6.0) * (dangVel1 + 2*dangVel2 + 2*dangVel3 + dangVel4)
        
        self.positionOffset = finalpos
        self.velocity = finalvelocity
        self.orientation = finalOrientation
        self.angularVelocity = finalangVel
        
        self.compute()

        if self.temperatureMap:
            self.history["maxTemperature"].append(float(np.max(self.temperatureMap)))
            self.history["averageTemperature"].append(float(np.mean(self.temperatureMap)))
        else:
            self.history["maxTemperature"].append(0.0)
            self.history["averageTemperature"].append(0.0)
        
        self.history["time"].append(float(Decimal(str(time)) + Decimal(str(dt))))
        self.history["position"].append(self.positionOffset.copy())
        self.history["velocity"].append(self.velocity.copy())
        self.history["force"].append(self.Force.copy())
        self.history["torque"].append(self.Torque.copy())
        self.history["orientation"].append(self.orientation.copy())
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
    def __init__(self, radius, a, resolution, thetaX, thetaY, material = Graphene, thickness=1e-6):
        self.radius = radius
        self.a = a
        self.resolution = resolution
        self.thickness = thickness
        self.material = material
        
        if abs(self.a) < 1e-15:
            area = np.pi * self.radius**2
        else:
            area = (np.pi / (6 * self.a**2)) * ((1 + 4 * self.a**2 * self.radius**2)**1.5 - 1)
        
        self.mass = self.material.density * area * self.thickness
        self.angularVelocity = np.array([0.0, 0.0, 0.0])
        self.angularAcceleration = np.array([0.0, 0.0, 0.0])
        
        dr = self.radius / self.resolution
        dphi = 2 * np.pi / self.resolution
        
        surfacePoints = []
        areaWeights = []
        
        for i in range(self.resolution):
            r = self.radius * (i + 0.5) / self.resolution
            
            for m in range(self.resolution):
                phi = 2 * np.pi * (m + 0.5) / self.resolution
                
                x = r * np.cos(phi)
                y = r * np.sin(phi)
                
                z = -self.a * (self.radius**2 - r**2)
                
                surfacePoints.append((x, y, z))
                
                dA = np.sqrt(1 + (2*self.a*x)**2 + (2*self.a*y)**2) * r * dr * dphi
                areaWeights.append(dA)
        
        surfacePoints = np.asarray(surfacePoints, dtype=float)
        areaWeights = np.asarray(areaWeights, dtype=float)
        
        totalWeight = np.sum(areaWeights)
        normalizedWeights = areaWeights / totalWeight
        
        self.centerOfMassLocal = np.sum(surfacePoints * normalizedWeights[:, None], axis=0)
        
        centeredPoints = surfacePoints - self.centerOfMassLocal
        identity = np.eye(3)
        
        self.inertiaTensor = np.zeros((3, 3))
        
        for points, weights in zip(centeredPoints, normalizedWeights):
            self.inertiaTensor += weights * (np.dot(points, points) * identity - np.outer(points, points))
        
        self.reflectivity = self.material.reflectivity
        self.absorptivity = self.material.absorptivity
        self.orientation = eul2quat(thetaX, thetaY) #w, x, y, z
        
        self.velocity = np.array([0.0, 0.0, 0.0])
        self.positionOffset = np.array([0.0, 0.0, 0.0])
        
        self.Force = None
        self.Torque = None
        
        self.xPoints=[]
        self.yPoints=[]
        self.zPoints=[]
        self.intensities = []
        self.temperatureMap = []
        
        self.history = {
            "time": [],
            "position": [],
            "velocity": [],
            "force": [],
            "torque": [],
            "orientation": [],
            "angularVelocity": [],
            "maxTemperature": [],
            "averageTemperature": []
        }
    
    def compute(self):
        self.xPoints=[]
        self.yPoints=[]
        self.zPoints=[]
        self.intensities = []
        
        self.temperatureMap = []
        Force = np.array([0.0, 0.0, 0.0])
        Torque = np.array([0.0, 0.0, 0.0])
        
        for i in range(self.resolution):
            for j in range(self.resolution):
                r = self.radius * (i + 0.5) / self.resolution
                phi = 2 * np.pi * (j + 0.5) / self.resolution
                
                x = r * np.cos(phi)
                y = r * np.sin(phi)
                
                z = -self.a * (self.radius**2 - r**2)
                
                surfacePosition = np.array([x, y, z])
                localPosition = surfacePosition - self.centerOfMassLocal
                
                localNormal = np.array([-2*self.a*x, -2*self.a*y, 1.0])
                normLength = np.linalg.norm(localNormal)
                
                if normLength > 1e-15:
                    localNormal = localNormal / normLength
                else:
                    localNormal = np.array([0.0, 0.0, 1.0])
                
                normal = rotateVector(localNormal, self.orientation)
                rotatedPos = rotateVector(localPosition, self.orientation)
                
                position = rotatedPos + self.positionOffset
                
                self.xPoints.append(position[0])
                self.yPoints.append(position[1])
                self.zPoints.append(position[2])
                
                I = gaussian(position[0], position[1], position[2])
                
                dr = self.radius / self.resolution
                dphi = 2 * np.pi / self.resolution
                
                dA = np.sqrt(1 + (2*self.a*x)**2 + (2*self.a*y)**2) * r * dr * dphi
                
                absorbedFlux = self.absorptivity * I
                temp = (absorbedFlux / (self.material.emissivity * sigma)) ** 0.25
                F = calcForce(I, normal, dA, self.reflectivity, self.absorptivity)
                
                tau = np.cross(rotatedPos, F)
                
                Force += F
                Torque += tau
                
                self.intensities.append(I)
                self.temperatureMap.append(temp)
                
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

    def derivatives(self, position, velocity, orientation, angularVelocity):
        oldPosition = self.positionOffset.copy()
        oldVelocity = self.velocity.copy()
        oldOrientation = self.orientation.copy()
        oldAngVel = self.angularVelocity.copy()
        
        self.positionOffset = position
        self.velocity = velocity
        self.orientation = orientation
        self.angularVelocity = angularVelocity
        
        self.compute()
        
        dpos = self.velocity.copy() #dx/dt = v
        dvel = self.Force / self.mass #dv/dt = a
        omegaQuat = np.array([0.0, angularVelocity[0], angularVelocity[1], angularVelocity[2]])
        qDot = 0.5 * quaternionMultiply(omegaQuat, orientation)
        # dangVel = np.linalg.solve(self.inertiaTensor, self.Torque) #d(omega)/dt = alpha # lmao "dang"
        
        R = quaternionToRotationMatrix(orientation)
        
        I_world = R @ self.inertiaTensor @ R.T
        angularMomentum = I_world @ angularVelocity
        
        dangVel = np.linalg.solve(I_world, self.Torque - np.cross(angularVelocity, angularMomentum))
        
        self.positionOffset = oldPosition
        self.velocity = oldVelocity
        self.orientation = oldOrientation
        self.angularVelocity = oldAngVel
        
        return dpos, dvel, qDot, dangVel
    
    def update(self, dt, time):
        pos = self.positionOffset.copy()
        velocity = self.velocity.copy()
        orientation = self.orientation.copy()
        angularVelocity = self.angularVelocity.copy()
        
        #k1
        dpos1, dvel1, dq1, dangVel1 = self.derivatives(pos, velocity, orientation, angularVelocity)
        
        #k2
        p2 = pos + 0.5 * dt * dpos1
        v2 = velocity + 0.5 * dt * dvel1
        q2 = orientation + 0.5 * dt * dq1
        q2 = normalizeQuaternion(q2)
        w2 = angularVelocity + 0.5 * dt * dangVel1
        
        dpos2, dvel2, dq2, dangVel2 = self.derivatives(p2, v2, q2, w2)
        
        #k3
        p3 = pos + 0.5 * dt * dpos2
        v3 = velocity + 0.5 * dt * dvel2
        q3 = orientation + 0.5 * dt * dq2
        q3 = normalizeQuaternion(q3)
        w3 = angularVelocity + 0.5 * dt * dangVel2
        
        dpos3, dvel3, dq3, dangVel3 = self.derivatives(p3, v3, q3, w3)
        
        #k4
        p4 = pos + dt * dpos3
        v4 = velocity + dt * dvel3
        q4 = orientation + dt * dq3
        q4 = normalizeQuaternion(q4)
        w4 = angularVelocity + dt * dangVel3
        
        dpos4, dvel4, dq4, dangVel4 = self.derivatives(p4, v4, q4, w4)
        
        #weighted average or something idk im just following a random RK4 tutorial i aint that smart
        finalpos = pos + (dt/ 6.0) * (dpos1 + 2*dpos2 + 2*dpos3 + dpos4)
        finalvelocity = velocity + (dt/ 6.0) * (dvel1 + 2*dvel2 + 2*dvel3 + dvel4)
        finalOrientation = orientation + (dt / 6.0) * (dq1 + 2*dq2 + 2*dq3 + dq4)
        finalOrientation = normalizeQuaternion(finalOrientation)
        finalangVel = angularVelocity + (dt/ 6.0) * (dangVel1 + 2*dangVel2 + 2*dangVel3 + dangVel4)
        
        self.positionOffset = finalpos
        self.velocity = finalvelocity
        self.orientation = finalOrientation
        self.angularVelocity = finalangVel

        if self.temperatureMap:
            self.history["maxTemperature"].append(float(np.max(self.temperatureMap)))
            self.history["averageTemperature"].append(float(np.mean(self.temperatureMap)))
        else:
            self.history["maxTemperature"].append(0.0)
            self.history["averageTemperature"].append(0.0)

        self.history["time"].append(float(Decimal(str(time)) + Decimal(str(dt))))
        self.history["position"].append(self.positionOffset.copy())
        self.history["velocity"].append(self.velocity.copy())
        self.history["force"].append(self.Force.copy())
        self.history["torque"].append(self.Torque.copy())
        self.history["orientation"].append(self.orientation.copy())
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
    time = Decimal('0.0')
    dt = timestep
    print("="*90)
    for step in range(steps):
        sail.compute()
        print(f"Position: {sail.positionOffset} \n Omega: {sail.angularVelocity} \n Torque: {sail.Torque} \n Force: {sail.Force} \n (Time: {time}) \n {'='*45}")
        sail.update(dt, time)
        time += Decimal(str(dt))
    sail.Visualize()
    plotForceHistory(sail)
    plotForceMagnitude(sail)
    # plotPositionHistory(sail)
    # plotVelocityHistory(sail)
    # plotAccelerationHistory(sail)
    # plotAngularVelocityHistory(sail)
    # plotTemperatureHistory(sail)
    # plotTemperatureHistogram(sail)
    # plotIntensityHistogram(sail)
    # plotKineticEnergyHistory(sail)
    # plotRotationalEnergyHistory(sail)
    # plotTotalEnergy(sail)
    # plotBeamOffset(sail)
    # plotEulerHistory(sail)

def e(resolution=100):
    sail = RectangleLightSail(10, 10, resolution, 0, 0, np.radians(0), np.radians(0))
    F, T = sail.compute()
    
    print(f"Max Temperature: {np.max(sail.temperatureMap)}")
    print(f"Average Temperature: {np.mean(sail.temperatureMap)}")

def f():
    sail = RectangleLightSail(10, 10, 500, 0, 0, np.radians(0), np.radians(0))
    F, T = sail.compute()
    
    print(f"Simulated Force: {F}")
    print(f"Simulated Torque: {T}")
    
    expectedForce = ((2 * sail.reflectivity + sail.absorptivity) /c) * (np.pi * I0 * w0**2/2)
    
    print(f"Expected Force: {expectedForce}")
    print(f"Force Error %: {np.linalg.norm(F[2] - expectedForce) / np.linalg.norm(expectedForce) * 100}%")

def bullshitAccel():
    sail1 = RectangleLightSail(10, 10, 100, 0, 0, np.radians(0), np.radians(0), thickness=1e-6)
    sail2 = RectangleLightSail(10, 10, 100, 0, 0, np.radians(0), np.radians(0), thickness=2e-6)
    
    sail1.compute()
    sail2.compute()
    
    a1 = sail1.Force / sail1.mass
    a2 = sail2.Force / sail2.mass
    
    print(f"Mass 1: {sail1.mass}\nMass 2: {sail2.mass}")
    print(f"Acceleration 1: {a1}\nAcceleration 2: {a2}")
    
    print(f"Acceleration Ratio: {np.linalg.norm(a1)/np.linalg.norm(a2)}")
    
    sail1.update(0.01, 0.0)
    sail1.Visualize()
    sail2.update(0.01, 0.0)
    sail2.Visualize()

def bullshitAngAccel():
    global beamX, beamY
    beamX=1
    sail1 = RectangleLightSail(10, 10, 100, 0, 0, np.radians(0), np.radians(0), thickness=1e-6)
    sail2 = RectangleLightSail(10, 10, 100, 0, 0, np.radians(0), np.radians(0), thickness=2e-6)

    sail1.compute()
    sail2.compute()

    alpha1 = np.linalg.solve(sail1.inertiaTensor, sail1.Torque)
    alpha2 = np.linalg.solve(sail2.inertiaTensor, sail2.Torque)

    print(f"Angular Acceleration 1: {alpha1}\nAngular Acceleration 2: {alpha2}")

    print(f"Angular Acceleration Ratio: {np.linalg.norm(alpha1)/np.linalg.norm(alpha2)}")
    
    sail1.update(0.01, 0.0)
    sail1.Visualize()
    sail2.update(0.01, 0.0)
    sail2.Visualize()

def lmaoTests():
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
    steps = 2
    # simLoop(dt, steps, nyxRect1)
    simLoop(dt, steps, nyxSphere1)
    # simLoop(dt, steps, nyxParaboloid1)
    # simLoop(dt, steps, nyxDisc1)
    # # for resolution in [10, 20, 50, 100, 200, 500]:
    # #     print(f"Running test with resolution {resolution}")
    # #     e(resolution=resolution)
    # # f()
    # bullshitAccel()
    # bullshitAngAccel()
    
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
    # print(f"Test Case 1: \n ====
    # ============== \n Total force is {TotalForce} \n Total torque is {TotalTorque} \n ==================")
    # nyxRect1.Visualize()

    # TotalForce, TotalTorque = nyxRect2.compute()
    # print(f"Test Case 2: \n ================== \n Total force is {TotalForce} \n Total torque is {TotalTorque} \n ==================")

    # TotalForce, TotalTorque = nyxRect3.compute()
    # print(f"Test Case 3: \n ================== \n Total force is {TotalForce} \n Total torque is {TotalTorque} \n ==================")

    # TotalForce, TotalTorque = nyxRect4.compute()
    # print(f"Test Case 4: \n ================== \n Total force is {TotalForce} \n Total torque is {TotalTorque} \n ==================")

if __name__ == '__main__':
    lmaoTests()