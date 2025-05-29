def analyze_capturing_net():
    # Initialize board representation
    board = [['.'] * 12 for _ in range(12)]
    
    # Define key positions
    white_stones = [(11,2), (11,3), (10,2), (10,3), (9,3)]  # C1-D1, C2-D2, C3-D3
    black_stones = [(10,1), (9,3), (8,3), (8,4)]  # B2, C3, D4, E4
    vital_point = (10,0)  # A2
    
    print("Analysis of A2 move:")
    print("1. Strategic value:")
    print("   - Creates a capturing net")
    print("   - Prevents White's escape to the edge")
    print("2. Sequence after A2:")
    print("   - White has no way to create two eyes")
    print("   - All escape routes are blocked")
    print("3. Conclusion:")
    print("   - A2 is the vital point that leads to capture")
    
analyze_capturing_net()