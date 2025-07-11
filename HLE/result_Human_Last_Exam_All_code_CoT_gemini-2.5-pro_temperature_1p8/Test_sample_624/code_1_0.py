import collections
import copy

def solve_rubiks_problem():
    """
    This function encapsulates the entire logic to solve the Rubik's Cube problem.
    It simulates the cube, applies the scramble, and searches for the shortest F2L solution.
    """

    # --- Step 1: Cube Representation and Move Logic ---
    
    def get_solved_cube():
        # U-D, L-R, F-B faces with W-Y, O-R, G-B colors
        return {
            'U': [['W'] * 3 for _ in range(3)], 'D': [['Y'] * 3 for _ in range(3)],
            'L': [['O'] * 3 for _ in range(3)], 'R': [['R'] * 3 for _ in range(3)],
            'F': [['G'] * 3 for _ in range(3)], 'B': [['B'] * 3 for _ in range(3)],
        }

    def apply_move(cube, move_str):
        new_cube = copy.deepcopy(cube)
        move_face = move_str[0]
        is_prime = "'" in move_str
        is_double = "2" in move_str
        
        rotations = 3 if is_prime else (2 if is_double else 1)

        for _ in range(rotations):
            if move_face == 'U':
                new_cube['U'] = [list(row) for row in zip(*new_cube['U'][::-1])]
                temp = new_cube['F'][0]
                new_cube['F'][0] = new_cube['R'][0]
                new_cube['R'][0] = new_cube['B'][0]
                new_cube['B'][0] = new_cube['L'][0]
                new_cube['L'][0] = temp
            elif move_face == 'D':
                new_cube['D'] = [list(row) for row in zip(*new_cube['D'][::-1])]
                temp = new_cube['F'][2]
                new_cube['F'][2] = new_cube['L'][2]
                new_cube['L'][2] = new_cube['B'][2]
                new_cube['B'][2] = new_cube['R'][2]
                new_cube['R'][2] = temp
            elif move_face == 'R':
                new_cube['R'] = [list(row) for row in zip(*new_cube['R'][::-1])]
                temp = [new_cube['U'][i][2] for i in range(3)]
                for i in range(3): new_cube['U'][i][2] = new_cube['F'][i][2]
                for i in range(3): new_cube['F'][i][2] = new_cube['D'][i][2]
                for i in range(3): new_cube['D'][i][2] = new_cube['B'][2 - i][0]
                for i in range(3): new_cube['B'][2 - i][0] = temp[i]
            elif move_face == 'L':
                new_cube['L'] = [list(row) for row in zip(*new_cube['L'][::-1])]
                temp = [new_cube['U'][i][0] for i in range(3)]
                for i in range(3): new_cube['U'][i][0] = new_cube['B'][2 - i][2]
                for i in range(3): new_cube['B'][2 - i][2] = new_cube['D'][i][0]
                for i in range(3): new_cube['D'][i][0] = new_cube['F'][i][0]
                for i in range(3): new_cube['F'][i][0] = temp[i]
            elif move_face == 'F':
                new_cube['F'] = [list(row) for row in zip(*new_cube['F'][::-1])]
                temp = new_cube['U'][2]
                for i in range(3): new_cube['U'][2][i] = new_cube['L'][2 - i][2]
                for i in range(3): new_cube['L'][i][2] = temp[i] if i == 1 else (new_cube['D'][0][0] if i == 0 else new_cube['D'][0][2])
                new_cube['L'][0][2], new_cube['L'][1][2], new_cube['L'][2][2] = new_cube['D'][0][2], new_cube['D'][0][1], new_cube['D'][0][0]
                for i in range(3): new_cube['D'][0][i] = new_cube['R'][2 - i][0]
                for i in range(3): new_cube['R'][i][0] = temp[i]
            elif move_face == 'B':
                new_cube['B'] = [list(row) for row in zip(*new_cube['B'][::-1])]
                temp = new_cube['U'][0]
                for i in range(3): new_cube['U'][0][i] = new_cube['R'][i][2]
                for i in range(3): new_cube['R'][i][2] = new_cube['D'][2][2-i]
                for i in range(3): new_cube['D'][2][i] = new_cube['L'][2 - i][0]
                for i in range(3): new_cube['L'][i][0] = temp[i]
        
        return new_cube
    
    def cube_to_tuple(cube):
        return tuple(tuple(map(tuple, cube[face])) for face in sorted(cube.keys()))

    # --- Step 2: Scramble the cube ---
    
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D"
    cube = get_solved_cube()
    for move in scramble.split():
        cube = apply_move(cube, move)

    # --- Step 3: F2L Pair Checking Logic (in reference to Johnny's orientation) ---
    
    def count_solved_f2l_pairs(c):
        solved = 0
        # J_FR -> WOG/OG -> Our UFL
        if c['U'][2][0]=='W' and c['F'][0][0]=='G' and c['L'][0][2]=='O' and c['F'][1][0]=='G' and c['L'][1][2]=='O': solved += 1
        # J_FL -> WOB/OB -> Our ULB
        if c['U'][2][2]=='W' and c['L'][0][0]=='O' and c['B'][0][2]=='B' and c['L'][1][0]=='O' and c['B'][1][2]=='B': solved += 1
        # J_BR -> WRG/RG -> Our URF
        if c['U'][0][2]=='W' and c['R'][0][0]=='R' and c['F'][0][2]=='G' and c['R'][1][0]=='R' and c['F'][1][2]=='G': solved += 1
        # J_BL -> WRB/RB -> Our URB
        if c['U'][0][0]=='W' and c['R'][0][2]=='R' and c['B'][0][0]=='B' and c['R'][1][2]=='R' and c['B'][1][0]=='B': solved += 1
        return solved
        
    # --- Step 4 & 5: BFS Search ---
    
    # Map Johnny's moves (Y-up, O-front) to our internal simulation's moves (W-up, G-front)
    move_map = {'U': 'D', 'D': 'U', 'F': 'L', 'L': 'B', 'B': 'R', 'R': 'F'}
    
    q = collections.deque([(cube, [])])
    visited = {cube_to_tuple(cube)}
    
    johnny_moves = [
        m + mod for m in "UDRLFB" for mod in ["", "'", "2"]
    ]

    while q:
        current_cube, path = q.popleft()
        
        if count_solved_f2l_pairs(current_cube) >= 2:
            sequence = " ".join(path)
            length = len(path)
            print("--- Solution Found ---")
            print(f"Minimum number of moves: {length}")
            # The prompt mentions outputting an equation. As this problem is about a move sequence,
            # this prints the solved sequence.
            print(f"Smallest sequence of moves: {sequence}")
            print(f"<<<{length}>>>")
            return

        if len(path) >= 8: # A reasonable limit for F2L searches
            continue
            
        for j_move in johnny_moves:
            our_move = move_map[j_move[0]] + j_move[1:]
            
            next_cube = apply_move(current_cube, our_move)
            next_cube_tuple = cube_to_tuple(next_cube)
            
            if next_cube_tuple not in visited:
                visited.add(next_cube_tuple)
                q.append((next_cube, path + [j_move]))
    
    print("No solution was found within the search limit.")

solve_rubiks_problem()