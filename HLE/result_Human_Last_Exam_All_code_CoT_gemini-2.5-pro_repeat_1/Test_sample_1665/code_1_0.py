def analyze_magnetic_scenarios():
    """
    Analyzes five scenarios to determine where a magnetic dipole's field
    will be strongest at the far end of a long cylinder.
    This is based on the physical properties of the materials.
    """

    print("Analyzing the magnetic field strength for five scenarios based on physical principles.")
    print("-" * 70)
    print("Key Principles:")
    print(" - Ferromagnetic Materials: Have high permeability and act as 'magnetic guides', concentrating field lines.")
    print(" - Ideal Superconductors: Exhibit the Meissner effect, expelling magnetic fields and acting as 'magnetic shields'.")
    print(" - Air (No Cylinder): Acts as a baseline where the field spreads out and weakens rapidly with distance.")
    print("-" * 70)

    scenarios = {
        1: {
            "description": "A solid ferromagnetic cylinder.",
            "analysis": "The high-permeability ferromagnetic material will capture the nearby magnetic field lines and guide them efficiently to the other end. This acts like a magnetic circuit yoke.",
            "result": "Strong magnetic field at the far end."
        },
        2: {
            "description": "A hollow tube made of an ideal superconducting material.",
            "analysis": "The superconducting material will expel the magnetic field (Meissner effect). Surface currents will shield the interior (the hollow core) and the region beyond the cylinder. The field is forced to go around the entire structure.",
            "result": "Very weak magnetic field at the far end."
        },
        3: {
            "description": "A ferromagnetic core surrounded by an ideal superconducting shell.",
            "analysis": "The outer superconducting shell is the first material the field encounters. It will expel the field, completely shielding the inner ferromagnetic core. The core never interacts with the field.",
            "result": "Very weak magnetic field at the far end."
        },
        4: {
            "description": "An ideal superconducting core surrounded by a ferromagnetic shell.",
            "analysis": "The outer ferromagnetic shell will capture and guide the magnetic field, just like in scenario 1. The field is concentrated within this shell. The inner superconducting core has little effect as the field is already being channeled around it.",
            "result": "Strong magnetic field at the far end."
        },
        5: {
            "description": "No cylinder, only air.",
            "analysis": "This is the baseline case. The dipole field spreads out in all directions and its strength decreases very rapidly with distance (as 1/r^3).",
            "result": "Very weak magnetic field at the far end."
        }
    }

    strongest_scenarios = []
    print("Scenario Analysis:")
    for i in range(1, 6):
        print(f"\nScenario {i}: {scenarios[i]['description']}")
        print(f"  - Analysis: {scenarios[i]['analysis']}")
        print(f"  - Result: {scenarios[i]['result']}")
        if "Strong" in scenarios[i]['result']:
            strongest_scenarios.append(i)
    
    print("-" * 70)
    print("Conclusion:")
    print("The magnetic field will be made stronger when a ferromagnetic material is positioned to capture and guide the field lines.")
    print("This occurs in scenarios 1 and 4.")
    
    # Final output of the identified scenario numbers
    print("\nThe final equation identifies the strongest field situations as:")
    equation_output = " and ".join(map(str, strongest_scenarios))
    print(f"Situations: {equation_output}")


if __name__ == "__main__":
    analyze_magnetic_scenarios()