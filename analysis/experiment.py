import csv
import json
import numpy as np

class Experiment:
    def __init__(self):
        self.results = []

    @staticmethod
    def _ensure_extension(filename, extension):
        filename = str(filename)
        if not filename.lower().endswith(extension.lower()):
            filename = f"{filename}{extension}"
        return filename
    
    def addResult(self, sail, geometry, parameters=None):
        result={
            "Geometry": geometry,
            "Parameters": parameters if parameters is not None else {},
            
            "FinalForce": sail.Force.tolist(),
            "FinalTorque": sail.Torque.tolist(),
            "FinalPosition": sail.positionOffset.tolist(),
            "FinalVelocity": sail.velocity.tolist(),
            "FinalAngularVelocity": sail.angularVelocity.tolist(),
            
            "Time": sail.history["time"],
            "PositionHistory": [p.tolist() for p in sail.history["position"]],
            "VelocityHistory": [v.tolist() for v in sail.history["velocity"]],
            "ForceHistory": [f.tolist() for f in sail.history["force"]],
            "TorqueHistory": [t.tolist() for t in sail.history["torque"]],
            "OrientationHistory": [q.tolist() for q in sail.history["orientation"]],
            "AngularVelocityHistory": [w.tolist() for w in sail.history["angularVelocity"]],
            
            "MaxTemperatureHistory": sail.history["maxTemperature"],
            "AverageTemperatureHistory": sail.history["averageTemperature"],
            
            "Mass": sail.mass,
            "Thickness": getattr(sail, "thickness", None),
            "MaxTemperature": max(sail.history["maxTemperature"]),
            "AverageTemperature": np.mean(sail.history["averageTemperature"])
        }
        self.results.append(result)
    
    def saveCSV(self, filename):
        if len(self.results) == 0:
            print("No result to Save.")
            return

        filename = self._ensure_extension(filename, ".csv")
        
        rows = []
        
        for result in self.results:
            row = {}
            for key, value in result["Parameters"].items():
                row[key] = value
                
            row["Geometry"] = result["Geometry"]
            row["Mass"] = result["Mass"]
            row["Thickness"] = result.get("Thickness", result["Parameters"].get("thickness", None))
            row["MaxTemperature"] = result["MaxTemperature"]
            row["AverageTemperature"] = result["AverageTemperature"]

            row["FinalPositionX"] = result["FinalPosition"][0]
            row["FinalPositionY"] = result["FinalPosition"][1]
            row["FinalPositionZ"] = result["FinalPosition"][2]
            
            row["FinalForceX"] = result["FinalForce"][0]
            row["FinalForceY"] = result["FinalForce"][1]
            row["FinalForceZ"] = result["FinalForce"][2]
            
            row["FinalTorqueX"] = result["FinalTorque"][0]
            row["FinalTorqueY"] = result["FinalTorque"][1]
            row["FinalTorqueZ"] = result["FinalTorque"][2]
            
            row["FinalVelocityX"] = result["FinalVelocity"][0]
            row["FinalVelocityY"] = result["FinalVelocity"][1]
            row["FinalVelocityZ"] = result["FinalVelocity"][2]
            
            row["FinalAngularVelocityX"] = result["FinalAngularVelocity"][0]
            row["FinalAngularVelocityY"] = result["FinalAngularVelocity"][1]
            row["FinalAngularVelocityZ"] = result["FinalAngularVelocity"][2]
            
            rows.append(row)
        
        fieldnames = sorted(set().union(*(row.keys() for row in rows)))
        
        with open(filename, "w", newline='') as file:
            writer = csv.DictWriter(file, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
    
    def saveJSON(self, filename):
        filename = self._ensure_extension(filename, ".json")
        with open(filename, "w") as file:
            json.dump(self.results, file, indent=4)
    
    def runSimulation(self, sail, geometry, duration, dt, parameters=None):
        time = 0.0
        
        while time < duration:
            sail.compute()
            
            sail.update(dt, time)
            time += dt
        
        sail.compute()
        
        self.addResult(sail, geometry, parameters)
    
    def runParameterSweep(self, geometryClass, geometryName, parameterList, duration, dt):
    
        for i, parameters in enumerate(parameterList, start=1): #9 13 21
            
            sail = geometryClass(**parameters)
            
            self.runSimulation(sail, geometryName, duration, dt, parameters)
            
            print(f"Finished {i}/{len(parameterList)} simulations.")
    
    def clearResults(self):
        self.results.clear()
    
    def numberOfResults(self):
        return len(self.results)
    
    def getResults(self, index):
        return self.results[index]
    
    def removeResult(self, index):
        del self.results[index]


