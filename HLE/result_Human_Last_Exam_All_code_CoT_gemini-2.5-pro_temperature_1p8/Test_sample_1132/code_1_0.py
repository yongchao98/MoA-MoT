import copy

class RubiksCube:
    """
    A class to represent and manipulate a 3x3 Rubik's Cube.
    The cube is represented by 6 faces, each a 3x3 matrix of colors.
    The orientation is fixed: F=Front, U=Up, R=Right, D=Down, L=Left, B=Back.
    """

    def __init__(self, initial_state):
        # The problem defines faces by color and orientation. We map them to standard notation.
        # F=White, U=Orange, R=Blue, L=Green, D=Red, B=Yellow
        self.faces = {
            'F': copy.deepcopy(initial_state['white']),
            'U': copy.deepcopy(initial_state['orange']),
            'R': copy.deepcopy(initial_state['blue']),
            'L': copy.deepcopy(initial_state['green']),
            'D': copy.deepcopy(initial_state['red']),
            'B': copy.deepcopy(initial_state['yellow']),
        }

    def _rotate_face_cw(self, face_key):
        """Rotates a face 90 degrees clockwise."""
        face = self.faces[face_key]
        self.faces[face_key] = [list(row) for row in zip(*face[::-1])]

    def _rotate_face_ccw(self, face_key):
        """Rotates a face 90 degrees counter-clockwise."""
        face = self.faces[face_key]
        self.faces[face_key] = [list(row)[::-1] for row in zip(*face)]

    def apply_move(self, move):
        """Applies a single move or a sequence of moves."""
        moves_map = {
            'R': self.R, 'U': self.U, 'F': self.F, 'L': self.L, 'D': self.D, 'B': self.B,
            "R'": self.R_p, "U'": self.U_p, "F'": self.F_p, "L'": self.L_p, "D'": self.D_p, "B'": self.B_p
        }
        if move in moves_map:
            moves_map[move]()
        else:
            raise ValueError(f"Unknown move: {move}")

    def F(self):
        self._rotate_face_cw('F')
        temp = [self.faces['U'][2][i] for i in range(3)]
        # L->U, D->L, R->D, U(temp)->R
        for i in range(3): self.faces['U'][2][i] = self.faces['L'][2-i][2]
        for i in range(3): self.faces['L'][i][2] = self.faces['D'][0][i]
        for i in range(3): self.faces['D'][0][i] = self.faces['R'][2-i][0]
        for i in range(3): self.faces['R'][i][0] = temp[i]

    def R(self):
        self._rotate_face_cw('R')
        temp = [self.faces['F'][i][2] for i in range(3)]
        # D->F, B->D, U->B, F(temp)->U
        for i in range(3): self.faces['F'][i][2] = self.faces['D'][i][2]
        for i in range(3): self.faces['D'][i][2] = self.faces['B'][2-i][0]
        for i in range(3): self.faces['B'][i][0] = self.faces['U'][i][2]
        for i in range(3): self.faces['U'][i][2] = temp[i]

    def U(self):
        self._rotate_face_cw('U')
        temp = self.faces['F'][0]
        # L->F, B->L, R->B, F(temp)->R
        self.faces['F'][0] = self.faces['L'][0]
        self.faces['L'][0] = self.faces['B'][0]
        self.faces['B'][0] = self.faces['R'][0]
        self.faces['R'][0] = temp

    def L(self):
        self._rotate_face_cw('L')
        temp = [self.faces['F'][i][0] for i in range(3)]
        # U->F, B->U, D->B, F(temp)->D
        for i in range(3): self.faces['F'][i][0] = self.faces['U'][i][0]
        for i in range(3): self.faces['U'][i][0] = self.faces['B'][2-i][2]
        for i in range(3): self.faces['B'][i][2] = self.faces['D'][i][0]
        for i in range(3): self.faces['D'][i][0] = temp[i]

    def D(self):
        self._rotate_face_cw('D')
        temp = self.faces['F'][2]
        # R->F, B->R, L->B, F(temp)->L
        self.faces['F'][2] = self.faces['R'][2]
        self.faces['R'][2] = self.faces['B'][2]
        self.faces['B'][2] = self.faces['L'][2]
        self.faces['L'][2] = temp
        
    def B(self):
        self._rotate_face_cw('B')
        temp = [self.faces['U'][0][i] for i in range(3)]
        # R->U, D->R, L->D, U(temp)->L
        for i in range(3): self.faces['U'][0][i] = self.faces['R'][i][2]
        for i in range(3): self.faces['R'][i][2] = self.faces['D'][2][2 - i]
        for i in range(3): self.faces['D'][2][i] = self.faces['L'][2 - i][0]
        for i in range(3): self.faces['L'][i][0] = temp[i]

    # Prime (counter-clockwise) moves are implemented as three clockwise moves
    def F_p(self): self.F(); self.F(); self.F()
    def R_p(self): self.R(); self.R(); self.R()
    def U_p(self): self.U(); self.U(); self.U()
    def L_p(self): self.L(); self.L(); self.L()
    def D_p(self): self.D(); self.D(); self.D()
    def B_p(self): self.B(); self.B(); self.B()

# Main execution
if __name__ == "__main__":
    # Define the initial jumbled state of the cube from the problem description
    initial_cube_state = {
        'white': [['R', 'Y', 'Y'], ['G', 'W', 'W'], ['G', 'O', 'O']],
        'orange': [['R', 'Y', 'W'], ['B', 'O', 'Y'], ['Y', 'R', 'O']],
        'blue': [['G', 'G', 'G'], ['R', 'B', 'B'], ['B', 'B', 'O']],
        'yellow': [['Y', 'W', 'R'], ['B', 'Y', 'O'], ['G', 'O', 'B']],
        'green': [['R', 'W', 'W'], ['R', 'G', 'O'], ['W', 'W', 'B']],
        'red': [['B', 'Y', 'Y'], ['R', 'R', 'G'], ['W', 'G', 'O']]
    }

    # Create a cube instance with the initial state
    cube = RubiksCube(initial_cube_state)

    # The 5-step algorithm to apply
    algorithm = ["R", "U", "F", "L'", "D"]

    # Apply each move in the algorithm
    for move in algorithm:
        cube.apply_move(move)

    # Get the final state of the White (Front) face
    final_white_face = cube.faces['F']

    # Print the final state of the white face
    for row in final_white_face:
        print(row)