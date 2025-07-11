import math

def solve_emperor_tomb_problem():
    """
    Calculates the optimal number of squares (N) and circles (M)
    to maximize the total number of engraved characters (K).
    """
    # Constants from the problem description
    SHEET_WIDTH = 140
    SHEET_HEIGHT = 110
    SQUARE_SIZE = 10
    CIRCLE_BBOX_SIZE = 40  # 20cm radius -> 40cm diameter
    CHARS_PER_SQUARE = 4
    
    # Calculation for characters per circle
    # Alphabet size = 7*7*9 = 441. Bagua symbols are base 8.
    # We need ceil(log8(441)) = 3 symbols per character.
    # Each circle holds 999 symbols. So, 999 / 3 = 333 characters per circle.
    CHARS_PER_CIRCLE = 333

    best_N = 0
    best_M = 0
    max_K = 0
    
    # This dictionary stores the results of the manual packing analysis.
    # For each number of circles (M), it gives the maximum number of squares (N)
    # that can fit in the remaining area of the 140x110cm sheet.
    # M: number of circles
    # N: max number of squares for that M
    packing_results = {
        0: 154, # M=0 (0x0 area), remaining 140x110 -> N = 14*11 = 154
        1: 138, # M=1 (40x40 area), remaining gives N = 138
        2: 122, # M=2 (80x40 area), remaining gives N = 122
        3: 106, # M=3 (120x40 area), remaining gives N = 106
        4: 90,  # M=4 (80x80 area), remaining gives N = 90
        5: 74,  # M=5 (complex L-shape), remaining gives N = 74
        6: 58   # M=6 (120x80 area), remaining gives N = 58
    }

    print("--- Evaluating Options ---")
    # Iterate through all possible numbers of circles (M)
    for M, N in packing_results.items():
        # Calculate the total number of characters for this combination
        K = (N * CHARS_PER_SQUARE) + (M * CHARS_PER_CIRCLE)
        print(f"Option: {M} circles (M) and {N} squares (N) -> Total Characters (K) = {K}")
        
        # Check if this is the new best result
        if K > max_K:
            max_K = K
            best_N = N
            best_M = M
            
    # --- Final Result ---
    chars_from_squares = best_N * CHARS_PER_SQUARE
    chars_from_circles = best_M * CHARS_PER_CIRCLE
    
    print("\n--- Optimal Solution Found ---")
    print(f"To maximize the engraved characters, we should produce {best_N} squares and {best_M} circles.")
    
    print("\n--- Calculation for Maximum Characters (K) ---")
    print(f"Number of squares (N): {best_N}")
    print(f"Number of circles (M): {best_M}")
    print(f"Characters per square: {CHARS_PER_SQUARE}")
    print(f"Characters per circle: {CHARS_PER_CIRCLE}")
    
    print("\nFinal Equation:")
    print(f"K = (N * Characters_per_square) + (M * Characters_per_circle)")
    print(f"K = ({best_N} * {CHARS_PER_SQUARE}) + ({best_M} * {CHARS_PER_CIRCLE})")
    print(f"K = {chars_from_squares} + {chars_from_circles}")
    print(f"K = {max_K}")

    print("\nFinal Answer (N:M:K):")
    print(f"{best_N}:{best_M}:{max_K}")

# Run the solver
solve_emperor_tomb_problem()