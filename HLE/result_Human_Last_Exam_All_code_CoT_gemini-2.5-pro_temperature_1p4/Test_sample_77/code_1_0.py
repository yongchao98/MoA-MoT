class SimulationComponent:
    """A base class for a component in our fluid simulation."""
    def __init__(self, name, component_type, is_present=True):
        self.name = name
        self.component_type = component_type
        self.is_present = is_present

    def __str__(self):
        return f"{self.name} ({self.component_type})"

def check_simulation_setup(components):
    """
    Checks if the essential components for a basic fluid simulation are present.
    In our case, we need 1 Domain, at least 1 Inflow, and at least 1 Obstacle.
    """
    component_counts = {"Domain": 0, "Inflow": 0, "Obstacle": 0}
    for comp in components:
        if comp.component_type in component_counts and comp.is_present:
            component_counts[comp.component_type] += 1

    # Check for the minimum required components
    domain_ok = component_counts["Domain"] == 1
    inflow_ok = component_counts["Inflow"] >= 1
    obstacle_ok = component_counts["Obstacle"] >= 1
    
    is_functional = domain_ok and inflow_ok and obstacle_ok

    print("--- Fluid Simulation Sanity Check ---")
    print(f"Checking for required components: Domain, Inflow, Obstacle...")
    print(f"Domain found: {component_counts['Domain']}")
    print(f"Inflow objects found: {component_counts['Inflow']}")
    print(f"Obstacle objects found: {component_counts['Obstacle']}")
    print("\nResult:")
    if is_functional:
        print("The simulation setup is functional. All required components are correctly configured.")
    else:
        print("The simulation setup is not functional. Missing one or more required components.")
    
    # Final equation showing the presence of each required component
    # We use 1 for present and 0 for absent in our "equation".
    # The sum must be 3 for a functional setup with these 3 components.
    print("\n--- Functional Equation ---")
    print(f"Presence(Domain) + Presence(Inflow) + Presence(Obstacle) = Total")
    domain_presence = 1 if domain_ok else 0
    inflow_presence = 1 if inflow_ok else 0
    obstacle_presence = 1 if obstacle_ok else 0
    total_score = domain_presence + inflow_presence + obstacle_presence
    print(f"{domain_presence} + {inflow_presence} + {obstacle_presence} = {total_score}")
    print(f"A score of 3 is required for this setup to be functional.")


# Representing the user's scene with our classes
domain = SimulationComponent("Simulation Box", "Domain")
fluid_emitter = SimulationComponent("Emitter Sphere", "Inflow")
collision_plane = SimulationComponent("Floor Obstacle", "Obstacle")

scene_components = [domain, fluid_emitter, collision_plane]

# Run the check
check_simulation_setup(scene_components)