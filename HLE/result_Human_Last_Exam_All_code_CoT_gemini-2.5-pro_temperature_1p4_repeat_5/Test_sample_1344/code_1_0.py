import textwrap

def suggest_catalyst():
    """
    This script outlines a proposed optimal catalyst system for both
    olefin polymerization and polyolefin hydrogenolysis.
    The choice is based on established principles from scientific literature,
    aiming for a robust, single-site system with dual functionality.
    """

    # --- Step 1: Explain the Core Challenge ---
    print("### Catalyst Design Plan ###")
    challenge_text = """
    The goal is to find a single-site catalyst that performs two opposing functions:
    1. Olefin Polymerization: Chain growth by adding monomers (e.g., ethylene, propylene).
       - Desired Conditions: Lower temperatures, high activity.
    2. Polyolefin Hydrogenolysis: Chain breaking into small alkanes (e.g., propane, butane).
       - Desired Conditions: Higher temperatures, presence of Hydrogen (H2).
    A successful catalyst must be 'schizophrenic' - stable and active under both sets of conditions.
    This suggests a supported catalyst for thermal stability.
    """
    print(textwrap.dedent(challenge_text))

    # --- Step 2: Select and Justify the Group IV Metal ---
    metal = {
        "element": "Zirconium (Zr)",
        "rationale": (
            "Zirconium is an excellent choice. It forms highly active polymerization catalysts and "
            "its hydride derivatives, especially when supported, are exceptionally effective for "
            "C-C bond cleavage. It strikes a good balance between the high reactivity of Titanium "
            "and the often higher cost/lower reactivity of Hafnium, making it a versatile and "
            "well-studied core for this dual purpose."
        )
    }
    print("--- 1. Optimal Metal Selection ---")
    print(f"Selected Metal: {metal['element']}")
    print("Rationale:")
    print(textwrap.fill(metal['rationale'], width=80))
    print("-" * 30)

    # --- Step 3: Select and Justify the Ligand ---
    ligand = {
        "name": "Constrained-Geometry Ligand (e.g., [η⁵-C₅Me₄(SiMe₂N-tBu)])",
        "rationale": (
            "A 'constrained-geometry' or 'half-sandwich' ligand provides the ideal electronic and "
            "steric environment. It has an open coordination sphere, which is crucial for both "
            "olefin insertion (polymerization) and accessing the polymer backbone C-C bonds (hydrogenolysis). "
            "This ligand class is famous for creating highly active, single-site polymerization catalysts."
        )
    }
    print("--- 2. Optimal Ligand Selection ---")
    print(f"Selected Ligand Type: {ligand['name']}")
    print("Rationale:")
    print(textwrap.fill(ligand['rationale'], width=80))
    print("-" * 30)

    # --- Step 4: Select and Justify the Support ---
    support = {
        "material": "Partially Dehydroxylated Silica (SiO₂)",
        "specification": "Calcined at 700°C (often denoted SiO₂-700)",
        "rationale": (
            "A solid support is critical for thermal stability during the high-temperature hydrogenolysis step. "
            "Grafting the Zirconium complex onto a silica surface creates isolated, truly single-site catalysts, "
            "preventing bimolecular decomposition. Using silica calcined at a high temperature (e.g., 700 °C) "
            "creates highly reactive siloxane bridge sites for grafting, leading to a more active and stable final catalyst."
        )
    }
    print("--- 3. Optimal Support Selection ---")
    print(f"Selected Support: {support['material']}")
    print(f"Support Specification: {support['specification']}")
    print("Rationale:")
    print(textwrap.fill(support['rationale'], width=80))
    print("-" * 30)

    # --- Step 5: Assemble and Print the Final Catalyst System ---
    print("\n### Final Proposed Catalyst System ###")
    print("Based on the analysis, the proposed optimal combination is:")

    final_catalyst_description = (
        f"A single-site {metal['element']} catalyst featuring a {ligand['name']} ligand, "
        f"grafted onto {support['material']} support that has been {support['specification']}."
    )
    print(textwrap.fill(final_catalyst_description, width=80))
    
    print("\nIn chemical notation, the active species on the surface would be represented as:")
    # This "equation" shows all the components as requested
    print(f"[{ligand['name']}]--{metal['element']}--H  (grafted onto)  {support['material']}-{support['specification']}")


if __name__ == '__main__':
    suggest_catalyst()