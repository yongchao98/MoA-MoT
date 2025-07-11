def find_shorter_lifetime_complexes():
    """
    Analyzes the stability of four Iridium complexes to determine which are
    expected to have shorter operational lifetimes in light-emitting devices.

    The key to solving this problem lies in understanding how fluorination
    of the phenylpyridine ligands affects the stability of the Iridium-Carbon bond.
    
    1.  General Rule: Fluorine is an electron-withdrawing group, which
        generally strengthens the Ir-C bond and increases complex stability.
        This makes Complex 4 (no fluorine) less stable than Complexes 1 and 2.

    2.  The Ortho-Effect: An important exception is when fluorine is at the
        ortho-position (position 2' of the phenyl ring). This configuration
        is known to create specific degradation pathways, making the complex
        significantly *less* stable. Complex 3 possesses this feature.

    Based on this, we can rank the stability (longest to shortest lifetime):
    Complex 1 > Complex 2 > Complex 4 > Complex 3.

    The complexes with shorter lifetimes are therefore Complex 3 and Complex 4.
    """
    
    # The complexes identified as having shorter lifetimes
    shorter_lifetime_complex_numbers = [3, 4]
    
    print("The complexes expected to show shorter lifetimes are:")
    # Using a loop to print each identified complex number clearly
    for number in sorted(shorter_lifetime_complex_numbers):
        print(f"Complex {number}")

# Execute the function to get the answer
find_shorter_lifetime_complexes()