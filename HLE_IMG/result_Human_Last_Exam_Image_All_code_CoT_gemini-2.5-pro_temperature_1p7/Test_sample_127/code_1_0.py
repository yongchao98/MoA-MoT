def analyze_ring_rotation():
    """
    This script analyzes the steric environment of each peripheral phenyl ring
    in the given Iridium complex to determine which can rotate freely.
    """
    
    # Ring attachments: {Ring Label: (Parent Ring, Position)}
    attachments = {
        'H': ('A', '4'),
        'J': ('C', '4'),
        'K': ('A', '2'),
        'L': ('C', '5'),
        'G': ('F', '6')
    }

    freely_rotating = []
    
    print("Step-by-step analysis of ring rotation:")

    # Analysis for H
    print("\nAnalyzing Ring H:")
    print("Ring H is at position 4 of ring A. This 'para' position is far from other bulky groups. Rotation is FREE.")
    freely_rotating.append('H')
    
    # Analysis for J
    print("\nAnalyzing Ring J:")
    print("Ring J is at position 4 of ring C. This 'para' position is also sterically unhindered. Rotation is FREE.")
    freely_rotating.append('J')
    
    # Analysis for K
    print("\nAnalyzing Ring K:")
    print("Ring K is at position 2 of ring A, which is 'ortho' to the bulky pyridine ring B. Rotation is RESTRICTED.")
    
    # Analysis for L
    print("\nAnalyzing Ring L:")
    print("Ring L is at position 5 of ring C, which is 'ortho' to the C-Ir bond and close to ligand E-F. Rotation is RESTRICTED.")

    # Analysis for G
    print("\nAnalyzing Ring G:")
    print("Ring G is at position 6 of ring F. It will clash with the adjacent ring E of the bipyridine ligand. Rotation is RESTRICTED.")
    
    # Final conclusion
    print("\n-------------------------")
    print("Final Conclusion:")
    print("The only rings that can rotate freely are those in sterically unhindered positions.")
    
    # The problem asks to output the labels in the final answer
    freely_rotating.sort() # To present in alphabetical order
    final_answer = ", ".join(freely_rotating)
    print(f"Freely rotating rings: {final_answer}")
    print("This corresponds to answer choice H.")

analyze_ring_rotation()
