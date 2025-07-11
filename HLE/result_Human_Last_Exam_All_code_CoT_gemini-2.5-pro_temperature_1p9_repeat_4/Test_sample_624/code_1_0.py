import sys
from collections import deque

def solve_f2l_pairs():
    """
    This script simulates a Rubik's cube to solve for the first two F2L pairs.
    1. It defines the cube state and all 18 possible moves (U/U'/U2, R/R'/R2, etc.).
    2. It applies the given scramble to a solved cube.
    3. It reorients the cube to have yellow on top and orange on front.
    4. It uses a Breadth-First Search (BFS) to find the shortest sequence of moves
       to solve any two of the four F2L pairs.
    5. Finally, it prints the length and the move sequence of the solution.
    """

    # --- Cube State and Move Definitions ---
    # The cube is represented as a tuple of 54 strings (colors).
    # Face order: U, R, F, D, L, B
    # Sticker indices:
    #        0  1  2
    #        3  4  5
    #        6  7  8
    # 36-44  18-26   9-17  45-53
    # (L)    (F)     (R)     (B)
    #       27-35
    #        (D)

    def move_U(state):
        s = list(state)
        # Rotate face
        s[0], s[1], s[2], s[3], s[5], s[6], s[7], s[8] = state[6], state[3], state[0], state[7], state[1], state[8], state[5], state[2]
        # Rotate sides
        s[18], s[19], s[20], s[9], s[10], s[11], s[45], s[46], s[47], s[36], s[37], s[38] = \
        state[36], state[37], state[38], state[18], state[19], state[20], state[9], state[10], state[11], state[45], state[46], state[47]
        return tuple(s)

    def move_R(state):
        s = list(state)
        s[9],s[10],s[11],s[12],s[14],s[15],s[16],s[17] = state[15],state[12],state[9],state[16],state[10],state[17],state[14],state[11]
        s[2],s[5],s[8],s[20],s[23],s[26],s[29],s[32],s[35],s[45],s[48],s[51] = \
        state[51],state[48],state[45],state[2],state[5],state[8],state[20],state[23],state[26],state[35],state[32],state[29]
        return tuple(s)

    def move_F(state):
        s = list(state)
        s[18],s[19],s[20],s[21],s[23],s[24],s[25],s[26] = state[24],state[21],state[18],state[25],state[19],state[26],state[23],state[20]
        s[6],s[7],s[8],s[17],s[14],s[11],s[27],s[28],s[29],s[36],s[39],s[42] = \
        state[42],state[39],state[36],state[6],state[7],state[8],state[17],state[14],state[11],state[29],state[28],state[27]
        return tuple(s)

    def move_D(state):
        s = list(state)
        s[27],s[28],s[29],s[30],s[32],s[33],s[34],s[35] = state[33],state[30],state[27],state[34],state[28],state[35],state[32],state[29]
        s[24],s[25],s[26],s[51],s[52],s[53],s[15],s[16],s[17],s[42],s[43],s[44] = \
        state[42],state[43],state[44],state[24],state[25],state[26],state[51],state[52],state[53],state[15],state[16],state[17]
        return tuple(s)

    def move_L(state):
        s = list(state)
        s[36],s[37],s[38],s[39],s[41],s[42],s[43],s[44] = state[42],state[39],state[36],state[43],state[37],state[44],state[41],state[38]
        s[0],s[3],s[6],s[47],s[50],s[53],s[33],s[30],s[27],s[18],s[21],s[24] = \
        state[18],state[21],state[24],state[0],state[3],state[6],state[47],state[50],state[53],state[33],state[30],state[27]
        return tuple(s)

    def move_B(state):
        s = list(state)
        s[45],s[46],s[47],s[48],s[50],s[51],s[52],s[53] = state[51],state[48],state[45],state[52],state[46],state[53],state[50],state[47]
        s[2],s[1],s[0],s[38],s[41],s[44],s[35],s[34],s[33],s[9],s[12],s[15] = \
        state[9],state[12],state[15],state[2],state[1],state[0],state[38],state[41],state[44],state[35],state[34],state[33]
        return tuple(s)

    moves = {
        'U': move_U, 'R': move_R, 'F': move_F,
        'D': move_D, 'L': move_L, 'B': move_B
    }

    # Generate prime and double moves
    for face in list(moves.keys()):
        moves[face + "'"] = lambda s, f=moves[face]: f(f(f(s)))
        moves[face + '2'] = lambda s, f=moves[f]: f(f(s))
    
    # --- Whole Cube Rotations ---
    def rotate_x(state): # R axis
        s = move_R(move_L(move_L(move_L(state))))
        # also rotate centers. M' move not implemented, so do manually
        m = list(s)
        m[4], m[22], m[31], m[49] = s[49], s[4], s[22], s[31]
        return tuple(m)

    def rotate_y(state): # U axis
        s = move_U(move_D(move_D(move_D(state))))
        m = list(s)
        m[22], m[13], m[49], m[40] = s[40], s[22], s[13], s[49]
        return tuple(m)

    # --- F2L Goal Check ---
    def check_f2l_solved(cube):
        solved_pairs = 0
        # Expected center colors after x2 y2 rotation
        U_c, R_c, F_c, D_c, L_c, B_c = 'Y', 'G', 'O', 'W', 'B', 'R'

        # Pair 1: Front-Right (Orange-Green)
        if (cube[35]==D_c and cube[26]==F_c and cube[17]==R_c and cube[25]==F_c and cube[14]==R_c):
            solved_pairs += 1
        # Pair 2: Front-Left (Orange-Blue)
        if (cube[33]==D_c and cube[24]==F_c and cube[42]==L_c and cube[23]==F_c and cube[38]==L_c):
            solved_pairs += 1
        # Pair 3: Back-Right (Red-Green)
        if (cube[29]==D_c and cube[53]==B_c and cube[15]==R_c and cube[52]==B_c and cube[12]==R_c):
            solved_pairs += 1
        # Pair 4: Back-Left (Red-Blue)
        if (cube[27]==D_c and cube[51]==B_c and cube[44]==L_c and cube[47]==B_c and cube[41]==L_c):
            solved_pairs += 1
        return solved_pairs

    # --- Main Logic ---

    # 1. Initial State
    scramble_str = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D"
    # Solved cube state: U(W)-R(R)-F(G)-D(Y)-L(O)-B(B)
    solved_state = tuple('W'*9 + 'R'*9 + 'G'*9 + 'Y'*9 + 'O'*9 + 'B'*9)
    
    # 2. Apply Scramble
    scrambled_cube = solved_state
    for move_str in scramble_str.split():
        scrambled_cube = moves[move_str](scrambled_cube)

    # 3. Reorient Cube: Yellow-top, Orange-front (x2 y2)
    reoriented_cube = rotate_x(rotate_x(scrambled_cube))
    reoriented_cube = rotate_y(rotate_y(reoriented_cube))

    # 4. BFS Search
    start_state = reoriented_cube
    
    if check_f2l_solved(start_state) >= 2:
        print("0 moves required. Two F2L pairs are already solved.")
        return

    q = deque([(start_state, [])]) # (state, path)
    visited = {start_state}
    
    move_names = list(moves.keys())
    
    # Map moves to their faces to prevent redundant sequences like R R' or F F2
    move_face = {name: name[0] for name in move_names}

    while q:
        current_state, path = q.popleft()
        
        last_move_face = move_face[path[-1]] if path else None

        for move in move_names:
            if last_move_face and move_face[move] == last_move_face:
                continue

            new_state = moves[move](current_state)

            if new_state not in visited:
                new_path = path + [move]
                if check_f2l_solved(new_state) >= 2:
                    print(f"The minimum length to solve two F2L pairs is: {len(new_path)}")
                    print("One possible sequence of moves is:")
                    # Outputting the equation with each number as requested
                    for i, m in enumerate(new_path):
                        print(m, end=" " if i < len(new_path) - 1 else "\n")
                    # Final answer format for parsing
                    print(f"\n<<<{len(new_path)}>>>")
                    return
                
                visited.add(new_state)
                q.append((new_state, new_path))

solve_f2l_pairs()