import collections

def solve_f2l_pairs():
    """
    This script simulates a Rubik's cube to find the shortest sequence of moves
    to solve two F2L pairs from a given scrambled state.
    """

    # --- Step 1: Cube Representation and Move Logic ---
    # The state of the cube is a tuple of 54 strings, representing sticker colors.
    # Order of faces: U, R, F, D, L, B
    # Colors for solved state: W, R, G, Y, O, B

    # Defines how stickers move for each clockwise face turn.
    # Each number is a sticker's index. (0-8:U, 9-17:R, 18-26:F, 27-35:D, 36-44:L, 45-53:B)
    # The cycles are defined as (from, to) sticker movements.
    cycles = {
        # Face, (face sticker cycle), (3 adjacent edge cycles)
        'U': [(0,2,8,6,1,5,7,3), (18,9,45,36), (19,10,46,37), (20,11,47,38)],
        'R': [(9,11,17,15,10,14,16,12), (8,26,35,47), (5,23,32,50), (2,20,29,53)],
        'F': [(18,20,26,24,19,23,25,21), (6,38,29,15), (7,41,28,12), (8,44,27,9)],
        'D': [(27,29,35,33,28,32,34,30), (24,17,51,44), (25,16,52,43), (26,15,53,42)],
        'L': [(36,38,44,42,37,41,43,39), (0,18,27,45), (3,21,30,48), (6,24,33,51)],
        'B': [(45,47,53,51,46,50,52,48), (2,36,33,11), (1,39,34,14), (0,42,35,17)],
    }

    def apply_move(state_tuple, move):
        """Applies a single move (e.g., 'R', "U'", 'F2') to a cube state."""
        state = list(state_tuple)
        face = move[0]
        is_prime = len(move) > 1 and move[1] == "'"

        move_cycles = cycles[face]
        
        # A prime move applies the cycles in reverse.
        if is_prime:
            for cycle in move_cycles:
                temp = state[cycle[-1]]
                for i in range(len(cycle) - 1, 0, -1):
                    state[cycle[i]] = state[cycle[i-1]]
                state[cycle[0]] = temp
        else: # Clockwise move
            for cycle in move_cycles:
                temp = state[cycle[0]]
                for i in range(len(cycle) - 1):
                    state[cycle[i]] = state[cycle[i+1]]
                state[cycle[-1]] = temp

        state_tuple = tuple(state)

        # A double move is just two clockwise moves.
        if len(move) > 1 and move[1] == '2':
            state_tuple = apply_move(state_tuple, face)

        return state_tuple

    # --- Step 2: Goal State Definition ---
    def count_solved_f2l_pairs(state):
        """
        Checks how many F2L pairs are solved for a white cross on the bottom (D) face.
        A pair is solved if its specific corner and edge pieces are in their
        home slots with the correct orientation.
        """
        solved_pairs = 0
        # Check F2L pair for Front-Right slot
        # Corner: White-Green-Red (WGR) at DFR. Edge: Green-Red (GR) at FR.
        # DFR stickers: D[8]=35(W), F[8]=26(G), R[2]=11(R).
        # FR stickers: F[5]=23(G), R[5]=14(R).
        if state[35]=='W' and state[26]=='G' and state[11]=='R' and \
           state[23]=='G' and state[14]=='R':
            solved_pairs += 1

        # Check F2L pair for Front-Left slot
        # Corner: White-Green-Orange (WGO) at DFL. Edge: Green-Orange (GO) at FL.
        # DFL stickers: D[6]=33(W), F[6]=24(G), L[2]=38(O).
        # FL stickers: F[3]=21(G), L[5]=41(O).
        if state[33]=='W' and state[24]=='G' and state[38]=='O' and \
           state[21]=='G' and state[41]=='O':
            solved_pairs += 1

        # Check F2L pair for Back-Right slot
        # Corner: White-Blue-Red (WBR) at DBR. Edge: Blue-Red (BR) at BR.
        # DBR stickers: D[2]=29(W), B[0]=45(B), R[8]=17(R).
        # BR stickers: B[1]=46(B), R[7]=16(R).
        if state[29]=='W' and state[45]=='B' and state[17]=='R' and \
           state[46]=='B' and state[16]=='R':
            solved_pairs += 1
            
        # Check F2L pair for Back-Left slot
        # Corner: White-Blue-Orange (WBO) at DBL. Edge: Blue-Orange (BO) at BL.
        # DBL stickers: D[0]=27(W), B[2]=47(B), L[8]=44(O).
        # BL stickers: B[3]=48(B), L[7]=43(O).
        if state[27]=='W' and state[47]=='B' and state[44]=='O' and \
           state[48]=='B' and state[43]=='O':
            solved_pairs += 1

        return solved_pairs

    # --- Step 3: Scramble Cube and Find Solution with BFS ---
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D"
    solved_state = tuple('W'*9 + 'R'*9 + 'G'*9 + 'Y'*9 + 'O'*9 + 'B'*9)

    # Apply the scramble to a solved cube
    scrambled_state = solved_state
    for move in scramble.split():
        scrambled_state = apply_move(scrambled_state, move)

    # After scrambling, Johnny flips the cube to put the white cross on bottom.
    # We don't need to rotate the cube state itself; the goal-checking function
    # `count_solved_f2l_pairs` already checks for the standard F2L case (white on bottom).
    start_state = scrambled_state

    # Check if the initial state already satisfies the condition
    if count_solved_f2l_pairs(start_state) >= 2:
        print("Two F2L pairs are already solved.")
        print("Sequence: (no moves needed)")
        return 0

    # BFS setup
    queue = collections.deque([(start_state, [])]) # (state, path_of_moves)
    visited = {start_state}
    all_moves = ["U", "D", "L", "R", "F", "B",
                 "U'", "D'", "L'", "R'", "F'", "B'",
                 "U2", "D2", "L2", "R2", "F2", "B2"]

    while queue:
        current_state, path = queue.popleft()

        for move in all_moves:
            next_state = apply_move(current_state, move)
            if next_state not in visited:
                if count_solved_f2l_pairs(next_state) >= 2:
                    final_path = path + [move]
                    path_str = " ".join(final_path)
                    print(f"The minimum sequence to solve two F2L pairs is:")
                    print(path_str)
                    print(f"The length of the sequence is: {len(final_path)}")
                    return len(final_path)

                visited.add(next_state)
                queue.append((next_state, path + [move]))

    print("No solution found within a reasonable depth.")
    return -1

if __name__ == '__main__':
    result_length = solve_f2l_pairs()
    print(f"<<<{result_length}>>>")