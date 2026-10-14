from analysis.experiment import Experiment
from SimMain import RectangleLightSail, SphereLightSail, ParaboloidLightSail

# & "C:\Users\dell\AppData\Local\Programs\Python\Python313\python.exe" "c:/Users/dell/Documents/Makeiums/Nyx/Simulation/smolTest1.py"

experiment = Experiment()

sail1 = RectangleLightSail(width=2, height=2, resolution=100, StartingX=0, StartingY=0, thetaX=0, thetaY=0)
sail2 = SphereLightSail(radius=2, resolution=100, thetaX=0, thetaY=0)
sail3 = ParaboloidLightSail(radius=2, a=3, resolution=100, thetaX=0, thetaY=0)

experiment.runSimulation(sail1, "Rect", 5, 0.1)
print(experiment.numberOfResults())
print(experiment.getResults(0))

experiment.saveCSV("SmolTest1-Rect.csv")
experiment.saveJSON("SmolTest1-Rect.json")
experiment.clearResults()

experiment.runSimulation(sail2, "Sphere", 5, 0.1)
print(experiment.numberOfResults())
print(experiment.getResults(0))

experiment.saveCSV("SmolTest1-Sphere.csv")
experiment.saveJSON("SmolTest1-Sphere.json")
experiment.clearResults()

experiment.runSimulation(sail3, "Paraboloid", 5, 0.1)
print(experiment.numberOfResults())
print(experiment.getResults(0))

experiment.saveCSV("SmolTest1-Parabola.csv")
experiment.saveJSON("SmolTest1-Parabola.json")
experiment.clearResults()