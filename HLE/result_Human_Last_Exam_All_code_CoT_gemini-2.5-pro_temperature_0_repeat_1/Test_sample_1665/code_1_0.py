def analyze_magnetic_scenarios():
    """
    Analyzes five scenarios involving a magnetic dipole and a cylinder
    to determine in which cases the magnetic field is strongest at the far end.
    """

    # Define the scenarios for clarity
    scenarios = {
        1: "A ferromagnetic cylinder.",
        2: "A hollow tube of ideal superconducting material.",
        3: "A ferromagnetic core with an ideal superconducting shell.",
        4: "A superconducting core with a ferromagnetic shell.",
        5: "No cylinder, only air (baseline)."
    }

    # Explanation of physical principles
    print("--- Magnetic Field Guiding Analysis ---")
    print("The strength of the magnetic field at the far end depends on the material's ability to guide magnetic flux.")
    print("\nKey Principles:")
    print(" - Ferromagnets (high permeability): Attract and concentrate magnetic field lines, acting as a 'magnetic pipe'.")
    print(" - Superconductors (zero permeability): Expel magnetic field lines (Meissner effect), acting as a perfect magnetic shield.")
    print(" - Air (baseline): Allows the field to spread out, causing it to weaken rapidly with distance.")

    # Analysis of each scenario
    print("\n--- Scenario Evaluation ---")
    analysis = {
        1: "STRONG GUIDING. The ferromagnetic material pulls in and channels the magnetic flux efficiently.",
        2: "STRONG GUIDING. The superconducting walls expel the field, forcing it to concentrate and travel down the hollow core.",
        3: "STRONGEST GUIDING. The ferromagnetic core gathers the flux, and the superconducting shell provides perfect shielding, preventing any leakage.",
        4: "MODERATE GUIDING. The ferromagnetic shell guides the flux, but the superconducting core expels it from the center, making the path less efficient.",
        5: "NO GUIDING. The field spreads out freely, resulting in the weakest field at the far end."
    }

    for i in range(1, 6):
        print(f"\nScenario {i}: {scenarios[i]}")
        print(f"  - Result: {analysis[i]}")

    # Final Conclusion
    print("\n--- Conclusion ---")
    print("Compared to the baseline case with only air (Scenario 5), any scenario with a cylinder that guides the magnetic field will result in a stronger field at the far end.")
    print("This is true for all scenarios involving ferromagnetic or superconducting materials.")
    
    stronger_scenarios = [1, 2, 3, 4]
    
    print("\nThe situations where the magnetic field will be stronger are:")
    # The final answer is the list of these scenarios.
    # We print each number as requested.
    print(f"Scenario: {stronger_scenarios[0]}")
    print(f"Scenario: {stronger_scenarios[1]}")
    print(f"Scenario: {stronger_scenarios[2]}")
    print(f"Scenario: {stronger_scenarios[3]}")

if __name__ == '__main__':
    analyze_magnetic_scenarios()
