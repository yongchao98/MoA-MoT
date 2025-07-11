def suggest_catalyst_system():
    """
    Suggests a candidate for a dual-function single-site catalyst for
    polyolefin polymerization and hydrogenolysis based on current
    catalysis research principles.
    """

    # --- 1. Define Catalyst Components ---
    metal_center = {
        "element": "Zirconium (Zr)",
        "group": "IV",
        "justification": "A classic metal for olefin polymerization. It can form highly reactive hydride species capable of sigma-bond metathesis, which is a potential pathway for C-C bond cleavage in hydrogenolysis."
    }

    ligand_system = {
        "type": "PNP Pincer Ligand",
        "example": "bis(2-diisopropylphosphino-4-methylphenyl)amide",
        "justification": "Provides exceptional thermal stability and creates a rigid, well-defined coordination environment. This prevents catalyst deactivation at high temperatures and allows precise tuning of the metal's electronic properties for dual reactivity."
    }

    support_material = {
        "material": "Mesoporous Silica (e.g., SBA-15)",
        "role": "Single-Site Immobilization",
        "justification": "A high-surface-area, inert support prevents the active metal centers from aggregating. The molecular catalyst is chemically grafted onto the silica surface, ensuring a uniform distribution of single active sites."
    }
    
    active_species = {
        "type": "Metal-Hydride (Zr-H)",
        "justification": "Formed in-situ with H2 gas, the Zr-H bond is the key active site for both initiating polymerization (via insertion) and breaking down polymer chains (via beta-alkyl elimination and hydrogenation)."
    }

    # --- 2. Print Catalyst Details and Justification ---
    print("--- Proposed Dual-Function Catalyst System ---")
    print(f"\n[Metal Center]: {metal_center['element']}")
    print(f"  Justification: {metal_center['justification']}")
    
    print(f"\n[Ligand System]: {ligand_system['type']}")
    print(f"  Justification: {ligand_system['justification']}")

    print(f"\n[Support Material]: {support_material['material']}")
    print(f"  Justification: {support_material['justification']}")

    print(f"\n[Proposed Active Species]: {active_species['type']}")
    print(f"  Justification: {active_species['justification']}")
    
    # --- 3. Output the Final "Equation" Components ---
    # The final equation represents the active catalyst species on its support.
    print("\n--- Final Catalyst Equation Components ---")
    
    component1 = f"({ligand_system['type']})"
    component2 = metal_center['element']
    component3 = "H" # Represents the active Hydride
    component4 = "@" # Represents 'supported on'
    component5 = support_material['material']
    
    print(f"Component 1 (Ligand): {component1}")
    print(f"Component 2 (Metal): {component2}")
    print(f"Component 3 (Active Group): {component3}")
    print(f"Component 4 (Anchor Symbol): {component4}")
    print(f"Component 5 (Support): {component5}")
    
    print("\nFull Catalyst Formula:")
    print(f"{component1}-{component2}-{component3} {component4} {component5}")


if __name__ == '__main__':
    suggest_catalyst_system()