try:
    from .materials import Material
except (ImportError, ModuleNotFoundError):
    from materials import Material

Graphene = Material(name='Graphene', reflectivity=0.95, absorptivity=0.05, emissivity=0.90, density=2200, maxTemperature=4500, colour='black')

Aluminium = Material(name='Aluminium', reflectivity=0.92, absorptivity=0.08, emissivity=0.06, density=2700, maxTemperature=933.47, colour='silver')

SiliconCarbide = Material(name='Silicon Carbide', reflectivity=0.88, absorptivity=0.12, emissivity=0.82, density=3210, maxTemperature=3000, colour='grey')

AllMaterials = [Graphene, Aluminium, SiliconCarbide]

if __name__ == "__main__":
    for material in AllMaterials:
        print(f"\nMaterial: {material.name}, Reflectivity: {material.reflectivity}, Absorptivity: {material.absorptivity}, Emissivity: {material.emissivity}, Density: {material.density} kg/m^3, Max Temperature: {material.maxTemperature} K, Colour: {material.colour} \n{"="*85}\n")