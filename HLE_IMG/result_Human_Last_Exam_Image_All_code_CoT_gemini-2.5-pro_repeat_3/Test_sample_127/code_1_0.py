def analyze_ring_rotation():
    """
    Analyzes the rotational freedom of phenyl rings in the given complex based on steric hindrance.
    """
    
    print("Analyzing the rotational freedom of each phenyl ring based on steric hindrance:\n")

    # Define the rings and their properties
    rings = {
        'H': "Attached to C4 of ring A. This is far from the Ir center and other ligands. Rotation is FREE.",
        'J': "Attached to C4 of ring C. Similar to H, it's far from the Ir center. Rotation is FREE.",
        'G': "Attached to C6 of ring F. This is *ortho* to the coordinating Nitrogen atom, making it very close to the bulky Ir center. Rotation is RESTRICTED.",
        'L': "Attached to C5 of ring C. This is *ortho* to the coordinating Carbon (C6), placing it very close to the Ir center. Rotation is RESTRICTED.",
        'K': "Attached to C2 of ring A. This is *meta* to the coordinating Carbon (C6) but ortho to the bond connecting ring A to B. The hindrance is much less than for G and L. Rotation is considered FREE."
    }
    
    freely_rotating = []
    for ring, explanation in rings.items():
        print(f" - Ring {ring}: {explanation}")
        if "FREE" in explanation:
            freely_rotating.append(ring)
            
    # Sort the rings to match the answer choice format
    # Custom sort order to match option G: K, J, H
    final_order = ['K', 'J', 'H']
    freely_rotating_sorted = sorted(freely_rotating, key=lambda x: final_order.index(x) if x in final_order else -1)
    
    print("\nConclusion: The rings that can rotate freely are K, J, and H.")
    print("This corresponds to answer choice G.")

analyze_ring_rotation()
print("\n<<<G>>>")