import collections

# U(W)-L(O)-F(G)-R(R)-B(B)-D(Y)
# W=0, O=1, G=2, R=3, B=4, Y=5
SOLVED_STATE = tuple(
    [0]*9 + [1]*9 + [2]*9 + [3]*9 + [4]*9 + [5]*9
)

# Sticker indices for face rotations
FACE_ROTATIONS = {
    'U': (0, 2, 8, 6, 1, 5, 7, 3), 'D': (45, 47, 53, 51, 46, 50, 52, 48),
    'F': (18, 20, 26, 24, 19, 23, 25, 21), 'B': (36, 38, 44, 42, 37, 41, 43, 39),
    'L': (9, 11, 17, 15, 10, 14, 16, 12), 'R': (27, 29, 35, 33, 28, 32, 34, 30),
}
# Sticker indices for side effects of moves
SIDE_ROTATIONS = {
    'U': ((18, 19, 20), (27, 28, 29), (36, 37, 38), (9, 10, 11)),
    'D': ((24, 25, 26), (15, 16, 17), (42, 43, 44), (33, 34, 35)),
    'F': ((6, 7, 8), (30, 33, 35), (47, 46, 45), (15, 12, 9)),
    'B': ((2, 1, 0), (11, 14, 17), (51, 52, 53), (32, 29, 27)),
    'L': ((0, 3, 6), (18, 21, 24), (45, 48, 51), (44, 41, 38)),
    'R': ((8, 5, 2), (42, 39, 36), (53, 50, 47), (26, 23, 20)),
}

def apply_move(state_tuple, move):
    state = list(state_tuple)
    face = move[0]
    turns = 1 if len(move) == 1 else (2 if move[1] == '2' else 3)
    
    for _ in range(turns):
        # Rotate face stickers
        f_map = FACE_ROTATIONS[face]
        prev_vals = [state[i] for i in f_map]
        for i in range(8):
            state[f_map[i]] = prev_vals[i-2]
            
        # Rotate side stickers
        s_map = SIDE_ROTATIONS[face]
        prev_vals = [state[i] for group in s_map for i in group]
        
        num_groups = len(s_map)
        group_size = len(s_map[0])
        for i in range(num_groups):
            for j in range(group_size):
                current_idx = i * group_size + j
                prev_idx = (current_idx - group_size) % len(prev_vals)
                state[s_map[i][j]] = prev_vals[prev_idx]

    return tuple(state)

def parse_scramble(scramble_str):
    moves = []
    i = 0
    while i < len(scramble_str):
        if i + 1 < len(scramble_str) and scramble_str[i+1] in ("'", "2"):
            moves.append(scramble_str[i:i+2])
            i += 2
        else:
            moves.append(scramble_str[i])
            i += 1
        # Skip space
        if i < len(scramble_str) and scramble_str[i] == ' ':
            i += 1
    return moves

# --- Goal check based on Johnny's orientation (Y-top, O-front) ---
# Our reference orientation is W-top(U), G-front(F)
# Johnny's D=Our U, J's F=Our L, J's R=Our B, J's L=Our F, J's B=Our R, J's U=Our D
# Colors: W=0, O=1, G=2, R=3, B=4
# Pair 1: J's DFR slot (W-O-B) is Our ULB slot
def check_wob(s):
    return (s[6]==0 and s[9+2]==1 and s[36+0]==4) and (s[9+5]==1 and s[36+1]==4)
# Pair 2: J's DFL slot (W-O-G) is Our ULF slot
def check_wog(s):
    return (s[8]==0 and s[9+0]==1 and s[18+0]==2) and (s[9+3]==1 and s[18+1]==2)
# Pair 3: J's DBR slot (W-R-B) is Our URB slot
def check_wrb(s):
    return (s[0]==0 and s[27+2]==3 and s[36+2]==4) and (s[R_edge_b_idx:=27+5]==3 and s[B_edge_r_idx:=36+3]==4)
# Pair 4: J's DBL slot (W-R-G) is Our URF slot
def check_wrg(s):
    return (s[2]==0 and s[27+0]==3 and s[18+2]==2) and (s[R_edge_f_idx:=27+3]==3 and s[F_edge_r_idx:=18+5]==2)

def count_solved_pairs(state):
    return sum([check_wob(state), check_wog(state), check_wrb(state), check_wrg(state)])

def main():
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D"
    
    current_state = SOLVED_STATE
    for move in parse_scramble(scramble):
        current_state = apply_move(current_state, move)
        
    # --- Verify white cross is solved in our reference frame ---
    # cross_ok = all([
    #     current_state[7]==0, current_state[9+1]==1, # U-L
    #     current_state[1]==0, current_state[18+7]==2, # U-F
    #     current_state[3]==0, current_state[27+1]==3, # U-R
    #     current_state[5]==0, current_state[36+7]==4  # U-B
    # ])
    # print(f"White cross solved after scramble: {cross_ok}") # It is.

    # BFS starts here
    queue = collections.deque([(current_state, [])])
    visited = {current_state}
    
    # Mapping Johnny's moves to our reference frame
    johnny_to_our_moves = {'U': 'D', 'D': 'U', 'F': 'L', 'L': 'F', 'B': 'R', 'R': 'B'}
    all_moves_johnny = [face + mod for face in "UDLFRB" for mod in ("", "'", "2")]

    if count_solved_pairs(current_state) >= 2:
        print("0")
        print("Two F2L pairs were already solved.")
        return

    while queue:
        state, path = queue.popleft()
        
        for move_j in all_moves_johnny:
            face_j = move_j[0]
            modifier = move_j[1:] if len(move_j) > 1 else ""
            
            # For U/D moves in Y-top, axis doesn't flip, so prime is prime.
            # For F/B, L/R they are on the equator of the reorientation, so direction flips
            # F(L) -> L(F), B(R)->R(B), F and L turn opposite directions. x2 reverses side turns
            if modifier == "'":
                 # x2 reverses L/R, y' reverses F/B and new L/R(old B/F)
                 # Overall, U/D/F/B/L/R Johnny maps to D'/U'/L/F/R/B our frame
                 mod_o = "" if face_j in "FBLR" else "'"
            elif modifier == "":
                 mod_o = "'" if face_j in "FBLR" else ""
            else: # mod == "2"
                 mod_o = "2"

            # Correction due to frame transformation complexities: U/D turn directions are preserved, others are inverted.
            if move_j[0] in "UD":
                move_o = johnny_to_our_moves[face_j] + modifier
            else: # FBLR
                mod_map = {"'":"", "": "'", "2":"2"}
                move_o = johnny_to_our_moves[face_j] + mod_map[modifier]
                
            new_state = apply_move(state, move_o)
            
            if new_state not in visited:
                new_path = path + [move_j]
                if count_solved_pairs(new_state) >= 2:
                    print("Solution found!")
                    print(f"Moves: {' '.join(new_path)}")
                    print(f"Length: {len(new_path)}")
                    # For the "equation format"
                    equation_str = " + ".join(["1"] * len(new_path))
                    print(f"Equation: {equation_str} = {len(new_path)}")
                    
                    final_length = len(new_path)
                    print(f"\n<<< {final_length} >>>")

                    return

                visited.add(new_state)
                queue.append((new_state, new_path))

if __name__ == '__main__':
    main()
