import sys

def analyze_covalency():
    """
    Analyzes and compares the covalency of two cerium compounds based on orbital overlap.
    """
    # Define the properties for the two cerium compounds based on the problem description.
    # We use numerical values to represent the key orbitals for the output.
    compound_1 = {
        "name": "CeF6(2-)",
        "metal_orbital": 4,  # Represents 4f
        "ligand_orbital": 2, # Represents 2p
        "overlap_strength": "greater" # Qualitative description from the problem
    }

    compound_2 = {
        "name": "CeCl6(2-)",
        "metal_orbital": 4,  # Represents 4f
        "ligand_orbital": 3, # Represents 3p
        "overlap_strength": "lesser"
    }

    # Step 1: State the fundamental principle.
    print("Principle: The degree of covalency in a bond is directly proportional to the extent of orbital overlap between the bonding atoms.")
    print("Greater orbital overlap leads to more significant electron sharing and thus stronger covalency.")
    print("-" * 50)

    # Step 2: Present the given information from the problem.
    # Note: We output each number in the final 'equation' or statement.
    print("Information given:")
    print(f"For {compound_1['name']}, the overlap is between the Ce {compound_1['metal_orbital']}f orbital and the F {compound_1['ligand_orbital']}p orbital.")
    print(f"For {compound_2['name']}, the overlap is between the Ce {compound_2['metal_orbital']}f orbital and the Cl {compound_2['ligand_orbital']}p orbital.")
    print(f"It is observed that the overlap in {compound_1['name']} is '{compound_1['overlap_strength']}'.")
    print("-" * 50)

    # Step 3: Draw the conclusion based on the principle.
    print("Conclusion:")
    if compound_1["overlap_strength"] == "greater":
        stronger_compound = compound_1["name"]
        weaker_compound = compound_2["name"]
    else:
        # This else case is for completeness but not triggered by the problem statement.
        stronger_compound = compound_2["name"]
        weaker_compound = compound_1["name"]

    print(f"Because the {compound_1['metal_orbital']}f-{compound_1['ligand_orbital']}p orbital overlap in {stronger_compound} is greater than the overlap in {weaker_compound},")
    print(f"{stronger_compound} will display stronger covalency.")

if __name__ == '__main__':
    analyze_covalency()