import collections

def get_moves():
    """
    Generates the 12 move functions for a Rubik's Cube.
    Each move function takes a state tuple and returns a new state tuple.
    The cube state is represented by a tuple of 54 integers, where each integer
    is a color ID for a sticker.
    Faces are ordered: U, L, F, R, B, D.
    """

    # --- Define sticker indices for each face ---
    # U: 0-8, L: 9-17, F: 18-26, R: 27-35, B: 36-44, D: 45-53
    
    # --- Define sticker permutations for the U (Up) move ---
    # These map the destination index to the source index.
    u_map = {
        # U face itself (clockwise)
        0: 6, 1: 3, 2: 0, 3: 7, 4: 4, 5: 1, 6: 8, 7: 5, 8: 2,
        # Side pieces
        9: 38, 10: 37, 11: 36, # B face top -> L face top
        18: 9, 19: 10, 20: 11, # L face top -> F face top
        27: 18, 28: 19, 29: 20, # F face top -> R face top
        36: 27, 37: 28, 38: 29, # R face top -> B face top
    }

    # --- Define whole-cube rotations to derive other moves from U ---
    # Rotation on X axis (like an R move)
    rot_x_map = {
        # U -> F
        18:0, 19:1, 20:2, 21:3, 22:4, 23:5, 24:6, 25:7, 26:8,
        # F -> D
        45:18, 46:19, 47:20, 48:21, 49:22, 50:23, 51:24, 52:25, 53:26,
        # D -> B (face reversed)
        44:53, 43:52, 42:51, 41:50, 40:49, 39:48, 38:47, 37:46, 36:45,
        # B -> U (face reversed)
        8:44, 7:43, 6:42, 5:41, 4:40, 3:39, 2:38, 1:37, 0:36,
        # L face (rotated)
        9:15, 12:16, 15:17, 10:12, 13:13, 16:14, 11:9, 14:10, 17:11,
        # R face (rotated)
        27:33, 30:34, 33:35, 28:30, 31:31, 34:32, 29:27, 32:28, 35:29
    }
    
    # Rotation on Y axis (like a U move)
    rot_y_map = {
        # L -> F -> R -> B -> L
        18:9, 19:10, 20:11, 21:12, 22:13, 23:14, 24:15, 25:16, 26:17,
        27:18, 28:19, 29:20, 30:21, 31:22, 32:23, 33:24, 34:25, 35:26,
        36:27, 37:28, 38:29, 39:30, 40:31, 41:32, 42:33, 43:34, 44:35,
        9:36, 10:37, 11:38, 12:39, 13:40, 14:41, 15:42, 16:43, 17:44,
        # U face (rotated)
        0:6, 1:3, 2:0, 3:7, 4:4, 5:1, 6:8, 7:5, 8:2,
        # D face (rotated)
        45:51, 46:48, 47:45, 48:52, 49:49, 50:46, 51:53, 52:50, 53:47
    }

    def _apply_map(state, amap):
        """Helper to apply a permutation map to a state."""
        new_state = list(state)
        for i in range(54):
            if i in amap:
                new_state[i] = state[amap[i]]
        return tuple(new_state)
    
    # --- Create move functions ---
    def U(state): return _apply_map(state, u_map)
    
    # Helper for creating other moves
    def _create_move(rotations, base_move_func):
        def move(state):
            s = state
            for r in rotations: s = r(s)
            s = base_move_func(s)
            # Apply inverse rotations
            for r in reversed(rotations): s = r(r(r(s)))
            return s
        return move

    # Rotations
    def rot_x(state): return _apply_map(state, rot_x_map)
    def rot_y(state): return _apply_map(state, rot_y_map)
    rot_z_map = rot_x_map # Just for creating the move, not a true Z rotation
    rot_z = _create_move([rot_y], rot_x) # Z = Y-rotated X
    
    # Base moves
    F = _create_move([rot_x], U)
    R = _create_move([rot_y], F)
    B = _create_move([rot_y], R)
    L = _create_move([rot_y], B)
    D = _create_move([rot_x, rot_x], U)

    # Prime moves (counter-clockwise)
    def U_prime(state): return U(U(U(state)))
    def F_prime(state): return F(F(F(state)))
    def R_prime(state): return R(R(R(state)))
    def B_prime(state): return B(B(B(state)))
    def L_prime(state): return L(L(L(state)))
    def D_prime(state): return D(D(D(state)))

    return [U, U_prime, L, L_prime, F, F_prime, R, R_prime, B, B_prime, D, D_prime]


def solve_rubiks_problem():
    """
    Calculates the number of permutations that solve the cube at move 4, 5, or 6.
    """
    # N[k] will store the number of ways to return to start in k steps.
    N = {}
    
    # Solved state: 9 stickers of color 0 (U), then 9 of 1 (L), etc.
    initial_state = tuple([i // 9 for i in range(54)])
    
    moves = get_moves()

    # `dists` stores: {state: number_of_paths_to_state}
    dists = {initial_state: 1}
    N[0] = 1

    # Perform a Breadth-First Search for 6 steps
    print("Calculating N(k) for k=1 to 6...")
    for k in range(1, 7):
        next_dists = collections.defaultdict(int)
        for state, count in dists.items():
            for move in moves:
                next_state = move(state)
                next_dists[next_state] += count
        
        dists = next_dists
        N[k] = dists.get(initial_state, 0)
        print(f"  N({k}) = {N[k]}")

    N4, N5, N6 = N[4], N[5], N[6]

    # The total number of successful permutations is given by the formula:
    # 132 * N(4) + 12 * N(5) + N(6)
    result = 132 * N4 + 12 * N5 + N6

    print("\nBased on the derived formula: |S4 U S5 U S6| = 132*N(4) + 12*N(5) + N(6)")
    print("The final calculation is:")
    print(f"132 * {N4} + 12 * {N5} + {N6} = {result}")

solve_rubiks_problem()
<<<24012>>>