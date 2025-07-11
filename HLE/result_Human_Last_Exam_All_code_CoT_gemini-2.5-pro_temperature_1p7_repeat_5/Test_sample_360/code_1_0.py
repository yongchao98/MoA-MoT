def count_isomers():
    """
    This function calculates the number of stereoisomers formed in the reaction between
    2,5-di(2-pyridyl)thiazolo[5,4-d]thiazole and cis-dichlorobis(bipyridine)ruthenium(II).
    
    The final product is a dinuclear complex of the form [(bpy)2Ru-(L)-Ru(bpy)2]^4+,
    where L is the symmetric bridging ligand with asymmetric chelating ends.
    """

    print("Step 1: Analyze the stereochemistry of the monomeric unit [Ru(bpy)2(Ptz)]^2+.")
    print("Each unit has two stereogenic elements:")
    print("  - Metal-centered chirality: Delta (Δ) or Lambda (Λ).")
    print("  - Ligand-orientation diastereomerism: A or B.")
    print("This gives 4 possible stereoisomeric units: ΔA, ΔB, ΛA, ΛB.")
    
    print("\nStep 2: Combine two units to form the dinuclear complex and count the isomers.")

    # Case 1: Homochiral metal centers (Δ,Δ or Λ,Λ)
    print("\n  Case a) Homochiral metal centers (both Δ or both Λ):")
    # [ΔA,ΔA] and [ΛA,ΛA] form one enantiomeric pair
    homo_pair_1 = 2
    # [ΔB,ΔB] and [ΛB,ΛB] form another enantiomeric pair
    homo_pair_2 = 2
    # [ΔA,ΔB] and [ΛA,ΛB] form a third enantiomeric pair
    homo_pair_3 = 2
    total_homochiral = homo_pair_1 + homo_pair_2 + homo_pair_3
    print(f"    - Combinations of same units ([ΔA,ΔA] / [ΛA,ΛA]): {homo_pair_1} isomers (1 enantiomeric pair)")
    print(f"    - Combinations of same units ([ΔB,ΔB] / [ΛB,ΛB]): {homo_pair_2} isomers (1 enantiomeric pair)")
    print(f"    - Combinations of different units ([ΔA,ΔB] / [ΛA,ΛB]): {homo_pair_3} isomers (1 enantiomeric pair)")
    print(f"    Subtotal for homochiral isomers = {total_homochiral}")
    
    # Case 2: Heterochiral metal centers (Δ,Λ)
    print("\n  Case b) Heterochiral metal centers (one Δ and one Λ):")
    # [ΔA,ΛA] is a meso compound (achiral)
    meso_1 = 1
    # [ΔB,ΛB] is another meso compound
    meso_2 = 1
    # [ΔA,ΛB] and [ΛA,ΔB] form an enantiomeric pair
    hetero_pair_1 = 2
    total_heterochiral = meso_1 + meso_2 + hetero_pair_1
    print(f"    - Meso compounds ([ΔA,ΛA] and [ΔB,ΛB]): {meso_1} + {meso_2} = {meso_1 + meso_2} isomers")
    print(f"    - Enantiomeric pair ([ΔA,ΛB] / [ΛA,ΔB]): {hetero_pair_1} isomers")
    print(f"    Subtotal for heterochiral isomers = {total_heterochiral}")
    
    # Final Calculation
    print("\nStep 3: Calculate the total number of isomers.")
    total_isomers = total_homochiral + total_heterochiral
    print(f"The total number of isomers is the sum of isomers from both cases.")
    print(f"Total = (Homochiral Isomers) + (Heterochiral Isomers)")
    print(f"Total = {total_homochiral} + {total_heterochiral} = {total_isomers}")
    
    print(f"\nFinal Answer: {total_isomers}")

count_isomers()