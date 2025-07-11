def calculate_odmr_resonances():
    """
    Calculates the number of ODMR resonances for a hypothetical NV center in a cubic lattice
    under an electric field parallel to a lattice edge.
    """
    # 1. Define the possible orientations of the NV center in a cubic lattice.
    # These are along the three mutually orthogonal Cartesian axes.
    nv_orientations = ['[100]', '[010]', '[001]']

    # 2. Define the E-field direction, parallel to one edge. Let's choose the z-axis ([001]).
    e_field_direction = '[001]'

    print(f"System: Hypothetical NV center in a cubic lattice.")
    print(f"Conditions: B=0, Electric field E is applied along the {e_field_direction} axis.\n")

    # 3. Use a set to store the unique resonance frequencies found.
    # We use symbolic names: 'D' for the zero-field splitting frequency and
    # 'delta' for the magnitude of the Stark shift.
    unique_frequencies = set()

    # 4. Analyze each orientation group.
    print("Analyzing the effect of the E-field on each NV orientation group:")
    print("----------------------------------------------------------------")

    for orientation in nv_orientations:
        # Case 1: NV axis is parallel to the E-field
        if orientation == e_field_direction:
            print(f"Group 1: NV axis parallel to E-field ({orientation})")
            print("  - The transverse component of the E-field is zero.")
            print("  - No Stark splitting of the m_s = +/-1 levels occurs.")
            print("  - This group contributes 1 resonance line at frequency: D\n")
            unique_frequencies.add("D")

        # Case 2: NV axis is perpendicular to the E-field
        else:
            print(f"Group 2/3: NV axis perpendicular to E-field ({orientation})")
            print("  - The E-field is purely transverse to the NV axis.")
            print("  - Stark effect splits the m_s = +/-1 levels.")
            print("  - This group contributes 2 resonance lines at frequencies: D - delta and D + delta\n")
            unique_frequencies.add("D - delta")
            unique_frequencies.add("D + delta")

    # 5. Summarize and count the unique frequencies.
    num_resonances = len(unique_frequencies)
    
    # Sort for display purposes. A simple sort doesn't handle the symbols well,
    # so we'll create a custom order.
    sorted_frequencies = []
    if "D - delta" in unique_frequencies: sorted_frequencies.append("D - delta")
    if "D" in unique_frequencies: sorted_frequencies.append("D")
    if "D + delta" in unique_frequencies: sorted_frequencies.append("D + delta")


    print("Summary:")
    print("----------------------------------------------------------------")
    print("The combined ODMR spectrum from all NV centers will show all unique frequencies.")
    print(f"The set of unique resonance frequencies is: {{{', '.join(sorted_frequencies)}}}")
    print("\nFinal Answer:")
    print(f"The total number of distinct resonances observed will be the number of unique frequencies.")
    print(f"Number of resonances = {num_resonances}")


if __name__ == '__main__':
    calculate_odmr_resonances()