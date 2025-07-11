import textwrap

def propose_catalyst_system():
    """
    Proposes a state-of-the-art bifunctional catalyst system for both
    olefin polymerization and polyolefin hydrogenolysis, and prints the rationale.
    """

    # --- Define the Catalyst Components ---
    group_iv_metal = "Hafnium (Hf)"
    ligand = "Constrained-geometry 'pincer' ligand (e.g., a PNP framework)"
    support = "Acidic solid support, specifically Sulfated Zirconia (S-ZrO2)"
    reactant = "Hydrogen gas (H2)"

    # --- Print the Explanation ---
    print("--- Proposal for an Optimal Bifunctional Catalyst System ---")
    print("\nThe challenge is to design a single-site catalyst that promotes two opposing reactions: ")
    print("1. Olefin Polymerization (building C-C bonds)")
    print("2. Polyolefin Hydrogenolysis (breaking C-C bonds in the presence of H2)")
    print("\nA promising approach is a tandem system where different parts of the catalyst perform different roles.")

    print("\n--- Proposed Catalyst Formulation ---")
    wrapper = textwrap.TextWrapper(initial_indent="   - ", subsequent_indent="     ")

    # 1. Group IV Metal
    print("\n1. Metal Center:")
    metal_desc = f"Metal: {group_iv_metal}. Hafnium is chosen from Group IV as it is known for producing high molecular weight polymers with high activity in single-site catalysis. It is also capable of activating H2 for the hydrogenolysis step."
    print(wrapper.fill(metal_desc))

    # 2. Ligand
    print("\n2. Ligand:")
    ligand_desc = f"Ligand: {ligand}. A robust pincer ligand provides thermal stability and creates a well-defined, single active site on the metal. Its specific geometry can leave a coordination site open for monomer/polymer binding while also facilitating the chemical steps of bond activation."
    print(wrapper.fill(ligand_desc))

    # 3. Support
    print("\n3. Support:")
    support_desc = f"Support: {support}. The support is not just an anchor but a co-catalyst. The strong Brønsted acid sites on Sulfated Zirconia are highly effective for cleaving the strong C-C bonds in the polyolefin backbone through a carbenium ion mechanism, which is essential for degradation into smaller alkanes."
    print(wrapper.fill(support_desc))
    
    # 4. Reactant
    print("\n4. Reactant for Degradation:")
    reactant_desc = f"Reactant: {reactant}. Hydrogen is required for the hydrogenolysis process. It caps the broken polymer fragments, resulting in the desired short, saturated alkanes."
    print(wrapper.fill(reactant_desc))

    print("\n--- Final Proposed Equation Components ---")
    print("The catalytic process doesn't have a simple, single equation. However, the key components of the proposed system are:")
    print(f"Group IV Metal: {group_iv_metal}")
    print(f"Ligand System: {ligand}")
    print(f"Active Support: {support}")
    print(f"Co-reactant for degradation: {reactant}")

# Execute the function to print the proposal
propose_catalyst_system()