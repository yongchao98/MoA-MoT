def identify_unstable_complexes():
    """
    This script analyzes four Iridium complexes to determine which are expected to
    have shorter operational lifetimes in light-emitting devices based on their
    molecular structure.
    """
    
    # Step 1: Define the complexes and their key structural features.
    # The key feature is the substitution pattern on the cyclometalating phenyl ring.
    complexes = {
        1: {'substituents': ['2\'-F (ortho)', '4\'-F (para)']},
        2: {'substituents': ['None']},
        3: {'substituents': ['2\'-F (ortho)']},
        4: {'substituents': ['4\'-F (para)']}
    }

    # Step 2: Explain the underlying chemical principle.
    print("Reasoning:")
    print("The lifetime of these iridium emitters is largely determined by their chemical stability.")
    print("A crucial factor for stability is the strength of the Iridium-Carbon (Ir-C) bond.")
    print("Electron-withdrawing groups like fluorine (F) on the phenyl ring weaken the Ir-C bond.")
    print("A weaker Ir-C bond leads to faster molecular degradation and a shorter device lifetime.")
    print("This destabilizing effect is particularly strong for fluorine at the 'ortho' position (next to the Ir-C bond).")
    print("-" * 50)
    
    # Step 3: Identify the complexes with shorter lifetimes.
    # We are looking for complexes with the most destabilizing feature, the ortho-fluorine.
    shorter_lifetime_complexes = []
    for complex_id, properties in complexes.items():
        # Check if any substituent string contains 'ortho'.
        has_ortho_fluorine = any('ortho' in sub for sub in properties['substituents'])
        if has_ortho_fluorine:
            shorter_lifetime_complexes.append(complex_id)
            
    shorter_lifetime_complexes.sort()
    
    # Step 4: Print the final conclusion and the identified numbers.
    print("Conclusion:")
    print("Complexes with ortho-fluorine substituents are expected to be the least stable and show shorter lifetimes.")
    print("The identified complexes are:")
    
    # As requested, output each number of the final answer.
    for num in shorter_lifetime_complexes:
        print(f"{num}")

# Execute the analysis
identify_unstable_complexes()