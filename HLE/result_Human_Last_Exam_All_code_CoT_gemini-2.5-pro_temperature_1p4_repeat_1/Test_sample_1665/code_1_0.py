def analyze_magnetic_scenarios():
    """
    Analyzes five scenarios to determine where a magnetic field will be strongest.
    """

    # 1. Define material properties
    materials = {
        "Ferromagnetic": "High permeability. Attracts and guides magnetic field lines, concentrating them.",
        "Superconductor": "Perfect diamagnet (zero permeability). Expels magnetic fields (Meissner effect).",
        "Air": "Low permeability (baseline). Allows magnetic fields to spread out."
    }

    # 2. Define scenarios and their analysis
    scenarios = {
        1: {
            "description": "Solid ferromagnetic cylinder.",
            "effect": "The entire cylinder acts as a magnetic guide. Field lines are strongly concentrated and channeled to the other end, resulting in a strong field.",
            "is_stronger": True
        },
        2: {
            "description": "Hollow ideal superconducting tube.",
            "effect": "The superconductor expels the magnetic field. The field is forced to travel around the outside of the tube, weakening it at the other end.",
            "is_stronger": False
        },
        3: {
            "description": "Ferromagnetic core, superconducting shell.",
            "effect": "The outer superconducting shell prevents the magnetic field from entering. The field is expelled, similar to case 2. The inner ferromagnetic core has no effect. The field will be weak.",
            "is_stronger": False
        },
        4: {
            "description": "Superconducting core, ferromagnetic shell.",
            "effect": "The outer ferromagnetic shell attracts and guides the magnetic field lines. The field is concentrated in the shell and channeled to the other end. This results in a strong field.",
            "is_stronger": True
        },
        5: {
            "description": "No cylinder (air only).",
            "effect": "This is the baseline case. The magnetic field spreads out, and its strength decreases significantly with distance.",
            "is_stronger": False # This is our reference, not stronger than itself.
        }
    }

    print("Analysis of Magnetic Field Strength in Different Scenarios\n")
    print("--- Key Material Properties ---")
    print(f"Ferromagnetic: {materials['Ferromagnetic']}")
    print(f"Superconductor: {materials['Superconductor']}\n")

    print("--- Scenario-by-Scenario Analysis ---")
    for i in sorted(scenarios.keys()):
        print(f"\nScenario {i}: {scenarios[i]['description']}")
        print(f"   Effect: {scenarios[i]['effect']}")

    # 3. Conclude and find the answer
    strong_scenarios = []
    for i, data in scenarios.items():
        if data["is_stronger"]:
            strong_scenarios.append(str(i))

    print("\n--- Conclusion ---")
    print("The situations where the magnetic field will be stronger are those where an outer ferromagnetic material guides the field.")
    print("These scenarios are:")
    # Final output of the numbers as requested
    for num in strong_scenarios:
        print(num)

# Execute the analysis
analyze_magnetic_scenarios()