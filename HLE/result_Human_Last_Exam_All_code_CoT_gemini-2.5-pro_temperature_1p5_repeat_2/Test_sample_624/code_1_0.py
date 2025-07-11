import collections
import copy

def solve_rubiks_f2l_problem():
    """
    Solves the Rubik's Cube F2L problem as described.
    This function encapsulates the entire logic: cube representation, moves,
    scrambling, reorientation, and the BFS search for the solution.
    """

    # Helper to rotate a 3x3 face matrix clockwise
    def rotate_face(face, times=1):
        new_face = copy.deepcopy(face)
        for _ in range(times):
            temp_face = copy.deepcopy(new_face)
            # Corners
            new_face[0][0], new_face[0][2], new_face[2][2], new_face[2][0] = \
                temp_face[2][0], temp_face[0][0], temp_face[0][2], temp_face[2][2]
            # Edges
            new_face[0][1], new_face[1][2], new_face[2][1], new_face[1][0] = \
                temp_face[1][0], temp_face[0][1], temp_face[1][2], temp_face[2][1]
        return new_face

    # Applies a move (e.g., 'R', "U'", "F2") to a given state
    def apply_move(state, move):
        s = copy.deepcopy(state)
        times = 1
        if "'" in move: times = 3
        elif "2" in move: times = 2
        
        move_char = move[0]
        
        for _ in range(times):
            # A full clockwise turn will be applied `times` times
            temp_s = copy.deepcopy(s)
            if move_char == 'U':
                s['U'] = rotate_face(temp_s['U'])
                s['L'][0], s['F'][0], s['R'][0], s['B'][0] = temp_s['F'][0], temp_s['R'][0], temp_s['B'][0], temp_s['L'][0]
            elif move_char == 'D':
                s['D'] = rotate_face(temp_s['D'])
                s['L'][2], s['B'][2], s['R'][2], s['F'][2] = temp_s['B'][2], temp_s['R'][2], temp_s['F'][2], temp_s['L'][2]
            elif move_char == 'R':
                s['R'] = rotate_face(temp_s['R'])
                for i in range(3):
                    s['U'][i][2] = temp_s['F'][i][2]
                    s['F'][i][2] = temp_s['D'][i][2]
                    s['D'][i][2] = temp_s['B'][2 - i][0]
                    s['B'][2 - i][0] = temp_s['U'][i][2]
            elif move_char == 'L':
                s['L'] = rotate_face(temp_s['L'])
                for i in range(3):
                    s['U'][i][0] = temp_s['B'][2 - i][2]
                    s['B'][2 - i][2] = temp_s['D'][i][0]
                    s['D'][i][0] = temp_s['F'][i][0]
                    s['F'][i][0] = temp_s['U'][i][0]
            elif move_char == 'F':
                s['F'] = rotate_face(temp_s['F'])
                for i in range(3):
                    s['U'][2][i] = temp_s['L'][2-i][2]
                    s['L'][i][2] = temp_s['D'][0][i]
                    s['D'][0][i] = temp_s['R'][2-i][0]
                    s['R'][i][0] = temp_s['U'][2][2-i]
                s['L'] = [list(col) for col in zip(*s['L'])] # Transpose to fix column update
                s['L'] = [list(col) for col in zip(*s['L'])]
            elif move_char == 'B':
                s['B'] = rotate_face(temp_s['B'])
                for i in range(3):
                    s['U'][0][i] = temp_s['R'][i][2]
                    s['R'][i][2] = temp_s['D'][2][2-i]
                    s['D'][2][i] = temp_s['L'][i][0]
                    s['L'][i][0] = temp_s['U'][0][2-i]
        return s

    # Convert state to a hashable tuple to store in the 'visited' set
    def state_to_tuple(state):
        return tuple(tuple(tuple(row) for row in state[face]) for face in sorted(state.keys()))

    # Function to check how many F2L pairs are solved for a Y-top, O-front orientation
    def count_solved_f2l_pairs(state):
        count = 0
        # Centers: U=Y, F=O, R=G, L=B, B=R, D=W
        # FR slot (Orange-Green)
        if (state['U'][2][2] == 'Y' and state['F'][0][2] == 'O' and state['R'][0][0] == 'G' and
            state['F'][1][2] == 'O' and state['R'][1][0] == 'G'):
            count += 1
        # FL slot (Orange-Blue)
        if (state['U'][2][0] == 'Y' and state['F'][0][0] == 'O' and state['L'][0][2] == 'B' and
            state['F'][1][0] == 'O' and state['L'][1][2] == 'B'):
            count += 1
        # BR slot (Red-Green)
        if (state['U'][0][2] == 'Y' and state['B'][0][0] == 'R' and state['R'][0][2] == 'G' and
            state['B'][1][0] == 'R' and state['R'][1][2] == 'G'):
            count += 1
        # BL slot (Red-Blue)
        if (state['U'][0][0] == 'Y' and state['B'][0][2] == 'R' and state['L'][0][0] == 'B' and
            state['B'][1][2] == 'R' and state['L'][1][0] == 'B'):
            count += 1
        return count

    # --- Main execution logic ---
    
    # 1. Define solved state and scramble sequence
    solved_state = {
        'U': [['W'] * 3 for _ in range(3)], 'L': [['O'] * 3 for _ in range(3)],
        'F': [['G'] * 3 for _ in range(3)], 'R': [['R'] * 3 for _ in range(3)],
        'B': [['B'] * 3 for _ in range(3)], 'D': [['Y'] * 3 for _ in range(3)],
    }
    scramble_seq = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D".split()
    all_moves = ['U', 'U\'', 'U2', 'D', 'D\'', 'D2', 'L', 'L\'', 'L2', 
                 'R', 'R\'', 'R2', 'F', 'F\'', 'F2', 'B', 'B\'', 'B2']

    # 2. Apply the scramble
    scrambled_state = solved_state
    for move in scramble_seq:
        scrambled_state = apply_move(scrambled_state, move)
        
    # 3. Reorient the cube (x2 y rotation)
    s = scrambled_state
    s_x2 = {
        'U': rotate_face(s['D'], 2), 'D': rotate_face(s['U'], 2),
        'F': rotate_face(s['B'], 2), 'B': rotate_face(s['F'], 2),
        'L': rotate_face(s['L'], 2), 'R': rotate_face(s['R'], 2)
    }
    initial_state = {
        'U': rotate_face(s_x2['U'], 1), 'D': rotate_face(s_x2['D'], 3),
        'F': s_x2['L'], 'L': s_x2['B'], 'B': s_x2['R'], 'R': s_x2['F']
    }

    # 4. Start BFS search
    if count_solved_f2l_pairs(initial_state) >= 2:
        print("Two F2L pairs are already solved.")
        print("0 = 0")
        return

    q = collections.deque([(initial_state, [])])
    visited = {state_to_tuple(initial_state)}

    while q:
        current_state, path = q.popleft()
        
        # To keep the search fast, we can prune if the path gets too long
        if len(path) >= 8:
            continue

        for move in all_moves:
            # Simple optimization: avoid redundant moves like U followed by U'
            if path and move[0] == path[-1][0]:
                continue
                
            next_state = apply_move(current_state, move)
            next_tuple = state_to_tuple(next_state)

            if next_tuple not in visited:
                new_path = path + [move]
                
                if count_solved_f2l_pairs(next_state) >= 2:
                    # Solution found
                    equation_parts = ["1"] * len(new_path)
                    print(f"Minimum number of moves: {len(new_path)}")
                    print(f"Solution sequence: {' '.join(new_path)}")
                    print(f"{' + '.join(equation_parts)} = {len(new_path)}")
                    return
                
                visited.add(next_tuple)
                q.append((next_state, new_path))
                
    print("No solution found within the search limit.")

solve_rubiks_f2l_problem()
<<<5>>>