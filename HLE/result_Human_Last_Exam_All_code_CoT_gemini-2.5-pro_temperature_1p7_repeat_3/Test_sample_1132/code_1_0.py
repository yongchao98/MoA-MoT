import copy

class RubiksCube:
    """
    A class to represent a Rubik's cube and its movements.
    The cube orientation is fixed: F(ront), B(ack), U(p), D(own), L(eft), R(ight).
    """

    def __init__(self, initial_state):
        self.faces = copy.deepcopy(initial_state)

    def _rotate_face_clockwise(self, face):
        return [[face[2-j][i] for j in range(3)] for i in range(3)]

    def _rotate_face_counter_clockwise(self, face):
        return [[face[j][2-i] for j in range(3)] for i in range(3)]

    def move_R(self):
        self.faces['R'] = self._rotate_face_clockwise(self.faces['R'])
        old_U, old_F, old_D, old_B = self.faces['U'], self.faces['F'], self.faces['D'], self.faces['B']
        for i in range(3):
            self.faces['U'][i][2], self.faces['F'][i][2], self.faces['D'][i][2], self.faces['B'][2-i][0] = \
            old_B[2-i][0], old_U[i][2], old_F[i][2], old_D[i][2]

    def move_U(self):
        self.faces['U'] = self._rotate_face_clockwise(self.faces['U'])
        old_F, old_R, old_B, old_L = self.faces['F'], self.faces['R'], self.faces['B'], self.faces['L']
        self.faces['F'][0], self.faces['R'][0], self.faces['B'][0], self.faces['L'][0] = \
        old_R[0], old_B[0], old_L[0], old_F[0]
        
    def move_F(self):
        self.faces['F'] = self._rotate_face_clockwise(self.faces['F'])
        old_U, old_R, old_D, old_L = self.faces['U'], self.faces['R'], self.faces['D'], self.faces['L']
        old_U_row = copy.deepcopy(old_U[2])
        for i in range(3):
            self.faces['U'][2][i] = old_L[2-i][2]
            self.faces['L'][2-i][2] = old_D[0][2-i]
            self.faces['D'][0][2-i] = old_R[i][0]
            self.faces['R'][i][0] = old_U_row[i]

    def move_L_prime(self):
        self.faces['L'] = self._rotate_face_counter_clockwise(self.faces['L'])
        old_U, old_F, old_D, old_B = self.faces['U'], self.faces['F'], self.faces['D'], self.faces['B']
        for i in range(3):
            self.faces['U'][i][0], self.faces['F'][i][0], self.faces['D'][i][0], self.faces['B'][2-i][2] = \
            old_F[i][0], old_D[i][0], old_B[2-i][2], old_U[i][0]
            
    def move_D(self):
        self.faces['D'] = self._rotate_face_clockwise(self.faces['D'])
        old_F, old_R, old_B, old_L = self.faces['F'], self.faces['R'], self.faces['B'], self.faces['L']
        self.faces['F'][2], self.faces['R'][2], self.faces['B'][2], self.faces['L'][2] = \
        old_L[2], old_F[2], old_R[2], old_B[2]
        
    def get_face(self, face_name):
        return self.faces[face_name]

# Define the initial jumbled state of the cube
initial_cube_state = {
    'F': [['R','Y','Y'],['G','W','W'],['G','O','O']], # White
    'U': [['R','Y','W'],['B','O','Y'],['Y','R','O']], # Orange
    'R': [['G','G','G'],['R','B','B'],['B','B','O']], # Blue
    'L': [['R','W','W'],['R','G','O'],['W','W','B']], # Green
    'D': [['B','Y','Y'],['R','R','G'],['W','G','O']], # Red
    'B': [['Y','W','R'],['B','Y','O'],['G','O','B']]  # Yellow
}

# Create a cube instance
my_cube = RubiksCube(initial_cube_state)

# Apply the 5-step algorithm
my_cube.move_R()
my_cube.move_U()
my_cube.move_F()
my_cube.move_L_prime()
my_cube.move_D()

# Get the final state of the white face (Front)
final_white_face = my_cube.get_face('F')

# Print the final face
print(final_white_face)

# To match the required output format, print each sticker value
print("\nFinal face values:")
for row in final_white_face:
    for sticker in row:
        print(sticker, end=' ')
    print()
