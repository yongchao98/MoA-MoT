def calculate_nmr_peaks():
    """
    Calculates the number of expected 1H NMR peaks for the given molecule based on its symmetry and structure.
    """
    
    # Molecule: 1,3,5-tri[((4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl)methyl]-2,4,6-trimethylbenzene

    # Step 1: Analyze the central benzene core.
    # The three methyl groups at positions 2, 4, and 6 are equivalent due to C3 symmetry.
    core_methyl_peaks = 1
    
    # Step 2: Analyze one of the three identical arms attached to the core.
    # The arm consists of a methylene linker and a camphor-indazole unit.
    
    # The two protons on the methylene (-CH2-) linker are diastereotopic.
    linker_methylene_peaks = 2
    
    # The camphor-indazole unit is chiral and rigid, leading to distinct signals for each non-equivalent proton.
    # - 1H on the pyrazole ring
    indazole_pyrazole_H_peaks = 1
    # - 1H at the bridgehead (C4)
    indazole_bridgehead_H_peaks = 1
    # - 2H at the C5 position (diastereotopic)
    indazole_C5_H_peaks = 2
    # - 2H at the C6 position (diastereotopic)
    indazole_C6_H_peaks = 2
    # - 3H from the methyl group at C7
    indazole_C7_methyl_peaks = 1
    # - 6H from the two diastereotopic methyl groups at C8
    indazole_C8_gem_dimethyl_peaks = 2
    
    # Sum the peaks from the indazole unit
    total_indazole_peaks = (indazole_pyrazole_H_peaks + 
                          indazole_bridgehead_H_peaks + 
                          indazole_C5_H_peaks + 
                          indazole_C6_H_peaks + 
                          indazole_C7_methyl_peaks + 
                          indazole_C8_gem_dimethyl_peaks)

    # Step 3: Sum all the peaks.
    total_peaks = core_methyl_peaks + linker_methylene_peaks + total_indazole_peaks

    # Print the detailed breakdown of the calculation.
    print("Analysis of expected 1H NMR peaks:")
    print("-" * 35)
    print(f"1. Central core methyls (3 equivalent -CH3): {core_methyl_peaks} peak")
    print(f"2. Methylene linker (-CH2-) protons (diastereotopic): {linker_methylene_peaks} peaks")
    print("3. Camphor-indazole moiety protons:")
    print(f"   - Pyrazole ring proton: {indazole_pyrazole_H_peaks} peak")
    print(f"   - Bridgehead proton (C4): {indazole_bridgehead_H_peaks} peak")
    print(f"   - C5 methylene protons: {indazole_C5_H_peaks} peaks")
    print(f"   - C6 methylene protons: {indazole_C6_H_peaks} peaks")
    print(f"   - C7 methyl protons: {indazole_C7_methyl_peaks} peak")
    print(f"   - C8 gem-dimethyl protons: {indazole_C8_gem_dimethyl_peaks} peaks")
    print("-" * 35)
    
    # The final equation requires showing each individual number being added.
    all_components = [
        core_methyl_peaks, linker_methylene_peaks, indazole_pyrazole_H_peaks,
        indazole_bridgehead_H_peaks, indazole_C5_H_peaks, indazole_C6_H_peaks,
        indazole_C7_methyl_peaks, indazole_C8_gem_dimethyl_peaks
    ]
    
    equation_str = " + ".join(map(str, all_components))
    print(f"Total number of peaks = {equation_str} = {total_peaks}")


# Run the calculation and print the result.
calculate_nmr_peaks()