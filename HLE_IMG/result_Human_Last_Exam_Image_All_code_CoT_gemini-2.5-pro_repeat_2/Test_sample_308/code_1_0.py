def solve_roller_puzzle():
    """
    Analyzes and matches roller configurations to their displacement plots.
    """
    
    # Step 1 & 2: Define features for configurations and plots
    # N_d: Number of lobes on the driver
    # wiggles: Number of wiggles in the plot
    # Other features: amplitude, symmetry, etc.
    
    analysis = {
        'A': {'wiggles': 2, 'notes': 'Asymmetric wiggles'},
        'B': {'wiggles': 4, 'notes': 'Large amplitude, seemingly identical wiggles'},
        'C': {'wiggles': 1, 'notes': 'One dominant wiggle'},
        'D': {'wiggles': 4, 'notes': '2+2 repeating pattern'},
        'E': {'wiggles': 3, 'notes': 'Identical wiggles'},
        'F': {'wiggles': 4, 'notes': 'Medium amplitude, seemingly identical wiggles'},
        'G': {'wiggles': 6, 'notes': '3x2 repeating pattern'},
        'H': {'wiggles': 4, 'notes': 'Small amplitude, seemingly identical wiggles'},
    }
    
    configs = {
        1: {'N_d': 4, 'notes': 'Asymmetric lobes (2+2), causes 2+2 wiggle pattern'},
        2: {'N_d': 6, 'notes': 'Interacts with 3-lobe driven, causes 3x2 wiggle pattern'},
        3: {'N_d': 2, 'notes': 'Asymmetric lobes'},
        4: {'N_d': 4, 'notes': 'Pointy lobes -> large amplitude'},
        5: {'N_d': 6, 'notes': 'Very shallow/smooth lobes -> small amplitude'},
        6: {'N_d': 4, 'notes': 'Rounded lobes -> medium amplitude'},
        7: {'N_d': 5, 'notes': 'One dominant lobe -> appears as 1 large wiggle'},
        8: {'N_d': 3, 'notes': 'Identical lobes'},
    }

    # Step 3: Perform the matching
    # This dictionary will store the final answer, mapping Plot letter to Config number.
    mapping = {}

    # Match by unique wiggle count
    mapping['A'] = 3  # Plot A (#w=2) -> Config 3 (N_d=2)
    mapping['E'] = 8  # Plot E (#w=3) -> Config 8 (N_d=3)
    mapping['G'] = 2  # Plot G (#w=6) -> Config 2 (N_d=6)

    # Match by unique pattern
    mapping['D'] = 1  # Plot D (2+2 pattern) -> Config 1 (asymmetric 4-lobes)
    
    # Match by exception
    mapping['C'] = 7  # Plot C (#w=1) -> Config 7 (N_d=5 with one dominant lobe)
    
    # Match remaining 4-wiggle plots (B, F, H) with remaining configs (4, 5, 6)
    # The counts mismatch: Config 5 has N_d=6, but must match a 4-wiggle plot.
    # We match by amplitude, assuming the smoothest config produces the lowest amplitude.
    # Amplitude order: B (large) > F (medium) > H (small)
    # Feature order: Config 4 (pointy) > Config 6 (rounded) > Config 5 (smoothest)
    mapping['B'] = 4
    mapping['F'] = 6
    mapping['H'] = 5

    # Step 4: Construct the final sequence
    print("Matching Logic:")
    print("Plot A (#w=2) matches Config 3 (Nd=2)")
    print("Plot B (large amp #w=4) matches Config 4 (pointy Nd=4)")
    print("Plot C (#w=1) matches Config 7 (dominant lobe Nd=5)")
    print("Plot D (2+2 pattern) matches Config 1 (asymmetric Nd=4)")
    print("Plot E (#w=3) matches Config 8 (Nd=3)")
    print("Plot F (med amp #w=4) matches Config 6 (rounded Nd=4)")
    print("Plot G (#w=6) matches Config 2 (Nd=6)")
    print("Plot H (small amp #w=4) matches Config 5 (smoothest Nd=6)")
    print("-" * 20)
    
    result_sequence = ""
    for plot_letter in sorted(mapping.keys()):
        config_number = mapping[plot_letter]
        result_sequence += str(config_number)

    print("Final sequence (Config number for Plots A-H):")
    print(result_sequence)

solve_roller_puzzle()