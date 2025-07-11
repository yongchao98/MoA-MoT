def get_point_group_info():
    """
    This function provides the symmetry point group for bis(2,5-dithiahexane)copper
    and displays its corresponding character table.
    """
    molecule_name = "bis(2,5-dithiahexane)copper"
    point_group = "D₂"

    # Explanation of the choice
    explanation = [
        f"Molecule: {molecule_name}",
        "1. Geometry: The complex consists of a central Copper(I) atom with two bidentate 2,5-dithiahexane ligands.",
        "   This results in a coordination number of 4, for which Cu(I) strongly prefers a tetrahedral geometry.",
        "2. Ideal Symmetry: An idealized tetrahedral M(A-A)₂ complex, where the ligands are symmetric, would have D₂d symmetry.",
        "3. Actual Symmetry: The ligand (CH₃-S-CH₂-CH₂-S-CH₃) is not planar due to the puckered ethylene bridge and the methyl groups.",
        "   This non-planarity eliminates the mirror planes (σ_d) present in the D₂d group.",
        "4. Conclusion: Removing the mirror planes from D₂d leaves the rotational subgroup D₂. This group contains the identity element (E) and three mutually perpendicular C₂ rotation axes.",
        f"\nTherefore, the point group for {molecule_name} is {point_group}.\n"
    ]

    # D2 Character Table
    character_table = {
        "Point Group": point_group,
        "Operations": ["E", "C₂(z)", "C₂(y)", "C₂(x)"],
        "A":          [1,   1,       1,       1],
        "B₁":         [1,   1,      -1,      -1],
        "B₂":         [1,  -1,       1,      -1],
        "B₃":         [1,  -1,      -1,       1],
    }

    # Printing the results
    print("\n".join(explanation))
    
    print(f"--- Character Table for {character_table['Point Group']} Point Group ---")
    
    # Print header
    header = f"{'':<4}" + "  ".join([f"{op:<7}" for op in character_table["Operations"]])
    print(header)
    print("-" * len(header))
    
    # Print rows
    for key in ["A", "B₁", "B₂", "B₃"]:
        row_str = f"{key:<4}"
        for val in character_table[key]:
            row_str += f"{val:<7d} "
        print(row_str)

# Execute the function to print the information
get_point_group_info()
