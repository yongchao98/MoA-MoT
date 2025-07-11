import copy

def solve():
    """
    Solves the Rubik's cube puzzle by simulating the given moves.
    """

    # 1. Represent the Cube and Face Mapping
    initial_cube_state = {
        'white': [['R', 'Y', 'Y'], ['G', 'W', 'W'], ['G', 'O', 'O']],  # Front
        'orange': [['R', 'Y', 'W'], ['B', 'O', 'Y'], ['Y', 'R', 'O']],  # Up
        'blue': [['G', 'G', 'G'], ['R', 'B', 'B'], ['B', 'B', 'O']],   # Right
        'yellow': [['Y', 'W', 'R'], ['B', 'Y', 'O'], ['G', 'O', 'B']],  # Back
        'green': [['R', 'W', 'W'], ['R', 'G', 'O'], ['W', 'W', 'B']],  # Left
        'red': [['B', 'Y', 'Y'], ['R', 'R', 'G'], ['W', 'G', 'O']]     # Down
    }

    # 2. Map standard notations to face colors
    faces = {
        'F': 'white', 'U': 'orange', 'R': 'blue', 'L': 'green', 'D': 'red', 'B': 'yellow'
    }

    def rotate_face_cw(matrix):
        """Rotates a 3x3 matrix 90 degrees clockwise."""
        return [list(row) for row in zip(*matrix[::-1])]

    # 3. Implement Rotation Functions
    def move_R(cube):
        """Performs a clockwise R (Right) move."""
        cube[faces['R']] = rotate_face_cw(cube[faces['R']])
        f, u, b, d = faces['F'], faces['U'], faces['B'], faces['D']
        temp_col = [cube[f][i][2] for i in range(3)]
        for i in range(3): cube[f][i][2] = cube[d][i][2]
        for i in range(3): cube[d][i][2] = cube[b][2-i][0]
        for i in range(3): cube[b][2-i][0] = cube[u][i][2]
        for i in range(3): cube[u][i][2] = temp_col[i]
        return cube

    def move_U(cube):
        """Performs a clockwise U (Up) move."""
        cube[faces['U']] = rotate_face_cw(cube[faces['U']])
        f, r, b, l = faces['F'], faces['R'], faces['B'], faces['L']
        temp_row = copy.deepcopy(cube[f][0])
        cube[f][0] = copy.deepcopy(cube[l][0])
        cube[l][0] = copy.deepcopy(cube[b][0])
        cube[b][0] = copy.deepcopy(cube[r][0])
        cube[r][0] = temp_row
        return cube

    def move_F(cube):
        """Performs a clockwise F (Front) move."""
        cube[faces['F']] = rotate_face_cw(cube[faces['F']])
        u, r, d, l = faces['U'], faces['R'], faces['D'], faces['L']
        temp_row = copy.deepcopy(cube[u][2])
        cube[u][2] = [cube[l][2-i][2] for i in range(3)]
        for i in range(3): cube[l][i][2] = cube[d][0][i]
        cube[d][0] = [cube[r][2-i][0] for i in range(3)]
        for i in range(3): cube[r][i][0] = temp_row[i]
        return cube

    def move_L(cube):
        """Performs a clockwise L (Left) move."""
        cube[faces['L']] = rotate_face_cw(cube[faces['L']])
        f, u, b, d = faces['F'], faces['U'], faces['B'], faces['D']
        temp_col = [cube[f][i][0] for i in range(3)]
        for i in range(3): cube[f][i][0] = cube[u][i][0]
        for i in range(3): cube[u][i][0] = cube[b][2-i][2]
        for i in range(3): cube[b][2-i][2] = cube[d][i][0]
        for i in range(3): cube[d][i][0] = temp_col[i]
        return cube

    def move_D(cube):
        """Performs a clockwise D (Down) move."""
        cube[faces['D']] = rotate_face_cw(cube[faces['D']])
        f, l, b, r = faces['F'], faces['L'], faces['B'], faces['R']
        temp_row = copy.deepcopy(cube[f][2])
        cube[f][2] = copy.deepcopy(cube[r][2])
        cube[r][2] = copy.deepcopy(cube[b][2])
        cube[b][2] = copy.deepcopy(cube[l][2])
        cube[l][2] = temp_row
        return cube

    # 4. Handle L'
    def move_L_prime(cube):
        """Performs a counter-clockwise L move (L' -> LLL)."""
        cube = move_L(cube)
        cube = move_L(cube)
        cube = move_L(cube)
        return cube

    # 5. Execute the Algorithm
    cube = copy.deepcopy(initial_cube_state)
    
    print(f"Initial White Face: {cube['white']}")

    cube = move_R(cube)
    print(f"After R move:       {cube['white']}")

    cube = move_U(cube)
    print(f"After U move:       {cube['white']}")

    cube = move_F(cube)
    print(f"After F move:       {cube['white']}")

    cube = move_L_prime(cube)
    print(f"After L' move:      {cube['white']}")

    cube = move_D(cube)
    print(f"After D move:       {cube['white']}")
    
    print("\nThe final state of the white face is:")
    print(cube['white'])

solve()
<<<A>>>