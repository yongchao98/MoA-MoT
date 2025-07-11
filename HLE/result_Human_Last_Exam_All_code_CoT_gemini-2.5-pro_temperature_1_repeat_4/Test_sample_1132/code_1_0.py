import copy

def rotate_cw(matrix):
    """Rotates a 3x3 matrix 90 degrees clockwise."""
    return [[matrix[2][0], matrix[1][0], matrix[0][0]],
            [matrix[2][1], matrix[1][1], matrix[0][1]],
            [matrix[2][2], matrix[1][2], matrix[0][2]]]

def rotate_ccw(matrix):
    """Rotates a 3x3 matrix 90 degrees counter-clockwise."""
    return [[matrix[0][2], matrix[1][2], matrix[2][2]],
            [matrix[0][1], matrix[1][1], matrix[2][1]],
            [matrix[0][0], matrix[1][0], matrix[2][0]]]

def apply_move(cube, move):
    """Applies a single move to the cube."""
    c = copy.deepcopy(cube)
    
    if move == 'R':
        c['R'] = rotate_cw(c['R'])
        temp_col = [c['U'][i][2] for i in range(3)]
        for i in range(3): c['U'][i][2] = c['F'][i][2]
        for i in range(3): c['F'][i][2] = c['D'][i][2]
        for i in range(3): c['D'][i][2] = c['B'][2-i][0]
        for i in range(3): c['B'][i][0] = temp_col[2-i]
    
    elif move == 'U':
        c['U'] = rotate_cw(c['U'])
        temp_row = c['F'][0]
        c['F'][0] = c['R'][0]
        c['R'][0] = c['B'][0]
        c['B'][0] = c['L'][0]
        c['L'][0] = temp_row
        
    elif move == 'F':
        c['F'] = rotate_cw(c['F'])
        temp_row = [c['U'][2][i] for i in range(3)]
        for i in range(3): c['U'][2][i] = c['L'][2-i][2]
        for i in range(3): c['L'][i][2] = c['D'][0][i]
        for i in range(3): c['D'][0][i] = c['R'][2-i][0]
        for i in range(3): c['R'][i][0] = temp_row[i]

    elif move == "L'":
        c['L'] = rotate_ccw(c['L'])
        temp_col = [c['F'][i][0] for i in range(3)]
        for i in range(3): c['F'][i][0] = c['D'][i][0]
        for i in range(3): c['D'][i][0] = c['B'][2-i][2]
        for i in range(3): c['B'][i][2] = c['U'][2-i][0]
        for i in range(3): c['U'][i][0] = temp_col[i]
        
    elif move == 'D':
        c['D'] = rotate_cw(c['D'])
        temp_row = c['F'][2]
        c['F'][2] = c['L'][2]
        c['L'][2] = c['B'][2]
        c['B'][2] = c['R'][2]
        c['R'][2] = temp_row
        
    return c

def main():
    # Initial jumbled state of the cube
    # F=Front(White), U=Up(Orange), R=Right(Blue), B=Back(Yellow), L=Left(Green), D=Down(Red)
    cube = {
        'F': [['R','Y','Y'],['G','W','W'],['G','O','O']],
        'U': [['R','Y','W'],['B','O','Y'],['Y','R','O']],
        'R': [['G','G','G'],['R','B','B'],['B','B','O']],
        'B': [['Y','W','R'],['B','Y','O'],['G','O','B']],
        'L': [['R','W','W'],['R','G','O'],['W','W','B']],
        'D': [['B','Y','Y'],['R','R','G'],['W','G','O']]
    }

    # Algorithm to execute
    algorithm = ['R', 'U', 'F', "L'", 'D']

    # Apply each move in the sequence
    for move in algorithm:
        cube = apply_move(cube, move)

    # Print the final state of the white face (Front)
    final_white_face = cube['F']
    print("Final White Face:")
    for row in final_white_face:
        print(row)

if __name__ == "__main__":
    main()