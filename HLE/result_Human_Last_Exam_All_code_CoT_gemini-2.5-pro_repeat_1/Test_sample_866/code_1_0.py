def solve_crease_pattern():
    """
    Determines the number of valid flat-foldable assignments for a given
    partially assigned crease pattern at a single vertex.
    """
    # The input from the problem description
    pattern_data = [60, 'M', 30, '?', 50, '?', 70, 'V', 150, '?']

    # Extract the crease assignments (M, V, ?) from the input list
    creases = [item for item in pattern_data if isinstance(item, str)]
    
    # Count the total number of creases
    N = len(creases)
    
    print("To determine the number of flat-foldable assignments, we analyze the conditions for single-vertex flat folding.")
    print("A necessary condition is Maekawa's Theorem, which relates the number of mountain (M) and valley (V) folds.")
    print("\nThe two main equations are:")
    print("1. M + V = N (where N is the total number of creases)")
    print("2. |M - V| = 2 (Maekawa's Theorem)")

    # Apply the rules to the given pattern
    print(f"\nFor the given pattern, the total number of creases is N = {N}.")
    print("Substituting this value, we must find integers M and V that satisfy the following system of equations:")
    
    # Output the specific equations for this problem, as requested
    print(f"M + V = {N}")
    print(f"|M - V| = 2")

    # Conclude based on the parity of N
    print("\nFor this system to have integer solutions for M and V, the total number of creases, N, must be an even number.")
    print(f"In this case, N = {N}, which is an odd number.")
    print("Therefore, it is impossible to find integer values for M and V that satisfy both equations.")
    print("This means that no matter how the unassigned folds ('?') are assigned, the pattern cannot be folded flat.")
    
    result = 0
    print(f"\nThe total number of valid assignments is: {result}")

solve_crease_pattern()