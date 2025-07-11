def count_nmr_peaks():
    """
    This function calculates the number of expected 1H NMR peaks for
    1,3,5-tri[((4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl)methyl]-2,4,6-trimethylbenzene
    based on chemical equivalence and molecular symmetry.
    """
    
    print("Analyzing the molecule to determine the number of non-equivalent proton signals...")
    print("-" * 70)
    
    # Part 1: Central Core (2,4,6-trimethylbenzene) and Linkers (-CH2-)
    # The molecule has C3 symmetry, making the three arms equivalent.
    
    # The three methyl groups on the central benzene ring (at C2, C4, C6) are equivalent.
    mesityl_methyl_peaks = 1
    
    # The three -CH2- linker groups are equivalent. However, the molecule is chiral
    # due to the (4S,7R) centers, making the two protons on each -CH2- group diastereotopic.
    # This results in two distinct signals for the CH2 protons.
    linker_ch2_peaks = 2
    
    core_and_linker_total = mesityl_methyl_peaks + linker_ch2_peaks
    
    # Part 2: One Substituent Arm
    # We analyze one ((4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl) group.
    # This is a rigid, chiral system, so many protons are non-equivalent.
    
    # Proton on the pyrazole ring (C3-H). Unique.
    pyrazole_proton_peaks = 1
    
    # Methyl groups on the bicyclic framework:
    # - One methyl at C7. Unique.
    c7_methyl_peaks = 1
    # - Two methyls at C8 (gem-dimethyl). Diastereotopic due to chirality.
    c8_gem_dimethyl_peaks = 2
    
    # Protons on the saturated rings of the framework:
    # - One proton at the C4 bridgehead. Unique.
    c4_bridgehead_proton_peaks = 1
    # - Two CH2 groups at C5 and C6. In the asymmetric framework, all four protons are non-equivalent.
    c5_c6_methylene_peaks = 4
    
    substituent_total = (pyrazole_proton_peaks + 
                         c7_methyl_peaks + 
                         c8_gem_dimethyl_peaks + 
                         c4_bridgehead_proton_peaks + 
                         c5_c6_methylene_peaks)

    # Part 3: Total Calculation
    total_peaks = core_and_linker_total + substituent_total
    
    print("The total number of peaks is the sum of peaks from each unique part:")
    print(f"1. Central mesityl methyl groups: {mesityl_methyl_peaks} peak")
    print(f"2. Linker -CH2- groups (diastereotopic): {linker_ch2_peaks} peaks")
    print("3. One substituent arm, which contains:")
    print(f"   - Pyrazole C3-H: {pyrazole_proton_peaks} peak")
    print(f"   - C7-methyl group: {c7_methyl_peaks} peak")
    print(f"   - C8 gem-dimethyl groups: {c8_gem_dimethyl_peaks} peaks")
    print(f"   - C4-bridgehead proton: {c4_bridgehead_proton_peaks} peak")
    print(f"   - C5/C6 methylene protons: {c5_c6_methylene_peaks} peaks")
    print(f"   -------------------------------------------")
    print(f"   Subtotal for one arm: {substituent_total} peaks")
    print("-" * 70)
    
    # Final Equation: Mesityl-Me + Linker-CH2 + Substituent Arm = Total
    print("Final Equation:")
    print(f"{mesityl_methyl_peaks} + {linker_ch2_peaks} + {substituent_total} = {total_peaks}")
    
    print(f"\nTotal number of expected 1H NMR peaks: {total_peaks}")

# Execute the function to get the answer.
count_nmr_peaks()