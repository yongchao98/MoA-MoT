import copy

class RubiksCube:
    def __init__(self, faces):
        # Deepcopy to prevent modifying the original input dictionary
        self.faces = copy.deepcopy(faces)

    def _rotate_face_cw(self, face_key):
        face = self.faces[face_key]
        self.faces[face_key] = [
            [face[2][0], face[1][0], face[0][0]],
            [face[2][1], face[1][1], face[0][1]],
            [face[2][2], face[1][2], face[0][2]]
        ]

    def _rotate_face_ccw(self, face_key):
        face = self.faces[face_key]
        self.faces[face_key] = [
            [face[0][2], face[1][2], face[2][2]],
            [face[0][1], face[1][1], face[2][1]],
            [face[0][0], face[1][0], face[2][0]]
        ]

    def move_R(self):
        self._rotate_face_cw('R')
        temp_col = [self.faces['U'][i][2] for i in range(3)]
        for i in range(3):
            self.faces['U'][i][2] = self.faces['F'][i][2]
            self.faces['F'][i][2] = self.faces['D'][i][2]
            self.faces['D'][i][2] = self.faces['B'][2-i][0]
            self.faces['B'][2-i][0] = temp_col[i]

    def move_U(self):
        self._rotate_face_cw('U')
        temp_row = self.faces['F'][0][:]
        self.faces['F'][0] = self.faces['R'][0][:]
        self.faces['R'][0] = self.faces['B'][0][:]
        self.faces['B'][0] = self.faces['L'][0][:]
        self.faces['L'][0] = temp_row

    def move_F(self):
        self._rotate_face_cw('F')
        temp_row = self.faces['U'][2][:]
        for i in range(3):
            self.faces['U'][2][i] = self.faces['L'][2-i][2]
            self.faces['L'][2-i][2] = self.faces['D'][0][2-i]
            self.faces['D'][0][2-i] = self.faces['R'][i][0]
            self.faces['R'][i][0] = temp_row[i]

    def move_L_prime(self):
        self._rotate_face_ccw('L')
        temp_col = [self.faces['F'][i][0] for i in range(3)]
        for i in range(3):
            self.faces['F'][i][0] = self.faces['D'][i][0]
            self.faces['D'][i][0] = self.faces['B'][2-i][2]
            self.faces['B'][2-i][2] = self.faces['U'][i][0]
            self.faces['U'][i][0] = temp_col[i]

    def move_D(self):
        self._rotate_face_cw('D')
        temp_row = self.faces['F'][2][:]
        self.faces['F'][2] = self.faces['L'][2][:]
        self.faces['L'][2] = self.faces['B'][2][:]
        self.faces['B'][2] = self.faces['R'][2][:]
        self.faces['R'][2] = temp_row

    def get_face(self, face_key):
        return self.faces[face_key]

def solve():
    # Initial orientation: F=White, U=Orange, R=Blue, B=Yellow, L=Green, D=Red
    initial_faces = {
        'F': [['R','Y','Y'],['G','W','W'],['G','O','O']], # White
        'U': [['R','Y','W'],['B','O','Y'],['Y','R','O']], # Orange
        'R': [['G','G','G'],['R','B','B'],['B','B','O']], # Blue
        'B': [['Y','W','R'],['B','Y','O'],['G','O','B']], # Yellow
        'L': [['R','W','W'],['R','G','O'],['W','W','B']], # Green
        'D': [['B','Y','Y'],['R','R','G'],['W','G','O']]  # Red
    }

    cube = RubiksCube(initial_faces)

    # 1. R
    cube.move_R()
    # 2. U
    cube.move_U()
    # 3. F
    cube.move_F()
    # 4. L'
    cube.move_L_prime()
    # 5. D
    cube.move_D()

    final_white_face = cube.get_face('F')
    print(final_white_face)

solve()