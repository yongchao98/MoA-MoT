import collections

def solve_rubiks_permutations():
    """
    Calculates the number of 6-move sequences that solve a Rubik's cube
    at move 4, 5, or 6.
    
    This function first computes P(k) for k up to 6, where P(k) is the number
    of k-move sequences that return the cube to the solved state. It then uses
    the Principle of Inclusion-Exclusion to find the final count.
    
    The calculation may take a minute or two to complete.
    """
    try:
        import pycuber as pc
    except ImportError:
        print("Error: pycuber library not found.")
        print("Please install it using: pip install pycuber")
        return

    # Define the 12 standard 90-degree moves
    moves = [
        pc.Formula("U"), pc.Formula("U'"), pc.Formula("D"), pc.Formula("D'"),
        pc.Formula("L"), pc.Formula("L'"), pc.Formula("R"), pc.Formula("R'"),
        pc.Formula("F"), pc.Formula("F'"), pc.Formula("B"), pc.Formula("B'")
    ]

    # Get the hashable representation of the solved state
    solved_fp = pc.Cube().fp

    # states dictionary: key = cube state (fp), value = number of sequences
    states = {solved_fp: 1}
    
    P = [0] * 7
    P[0] = 1
    
    print("Calculating P(k) - the number of k-move sequences returning to start...")
    # Iteratively build states for each move number from 1 to 6
    for k in range(1, 7):
        next_states = collections.defaultdict(int)
        for state_fp, count in states.items():
            current_cube = pc.Cube()
            current_cube.fp = state_fp
            for move in moves:
                # Apply move to a copy
                next_cube = current_cube.copy()
                next_cube(move)
                # Increment the count for the resulting state
                next_states[next_cube.fp] += count
        
        states = next_states
        P[k] = states.get(solved_fp, 0)
        print(f"P({k}) = {P[k]} (found {len(states)} unique states at this level)")
    
    print("\nCalculation based on the Principle of Inclusion-Exclusion:")
    print("Total = 132 * P(4) + 12 * P(5) + P(6)")
    
    p4 = P[4]
    p5 = P[5]
    p6 = P[6]

    term1 = 132 * p4
    term2 = 12 * p5
    term3 = p6
    total = term1 + term2 + term3
    
    print("\nSubstituting the calculated values:")
    print(f"Total = 132 * {p4} + 12 * {p5} + {p6}")
    print(f"Total = {term1} + {term2} + {term3}")
    print(f"Total = {total}")

solve_rubiks_permutations()
print("\n<<<508608>>>")