import copy

class RubiksCube:
    """
    A class to represent and manipulate a Rubik's Cube.
    """
    def __init__(self, initial_state):
        """
        Initializes the cube with a given state.
        The state is a dictionary mapping face names (F, U, R, L, B, D)
        to 3x3 matrices of colors.
        """
        self.faces = copy.deepcopy(initial_state)

    def _rotate_face_cw(self, face):
        """Rotates a face 90 degrees clockwise."""
        self.faces[face] = [list(row) for row in zip(*self.faces[face][::-1])]

    def _rotate_face_ccw(self, face):
        """Rotates a face 90 degrees counter-clockwise."""
        self.faces[face] = [list(row) for row in list(zip(*self.faces[face]))[::-1]]

    def R(self):
        """Performs a clockwise R move."""
        self._rotate_face_cw('R')
        temp = [self.faces['F'][i][2] for i in range(3)]
        for i in range(3): self.faces['F'][i][2] = self.faces['U'][i][2]
        for i in range(3): self.faces['U'][i][2] = self.faces['B'][2 - i][0]
        for i in range(3): self.faces['B'][2 - i][0] = self.faces['D'][i][2]
        for i in range(3): self.faces['D'][i][2] = temp[i]

    def U(self):
        """Performs a clockwise U move."""
        self._rotate_face_cw('U')
        temp = self.faces['F'][0]
        self.faces['F'][0] = self.faces['R'][0]
        self.faces['R'][0] = self.faces['B'][0]
        self.faces['B'][0] = self.faces['L'][0]
        self.faces['L'][0] = temp

    def F(self):
        """Performs a clockwise F move."""
        self._rotate_face_cw('F')
        temp = [row for row in self.faces['U'][2]]
        for i in range(3): self.faces['U'][2][i] = self.faces['L'][2 - i][2]
        for i in range(3): self.faces['L'][i][2] = self.faces['D'][0][i]
        for i in range(3): self.faces['D'][0][i] = self.faces['R'][2 - i][0]
        for i in range(3): self.faces['R'][i][0] = temp[i]

    def L_prime(self):
        """Performs a counter-clockwise L' move."""
        self._rotate_face_ccw('L')
        temp = [self.faces['F'][i][0] for i in range(3)]
        for i in range(3): self.faces['F'][i][0] = self.faces['U'][i][0]
        for i in range(3): self.faces['U'][i][0] = self.faces['B'][2 - i][2]
        for i in range(3): self.faces['B'][2 - i][2] = self.faces['D'][i][0]
        for i in range(3): self.faces['D'][i][0] = temp[i]

    def D(self):
        """Performs a clockwise D move."""
        self._rotate_face_cw('D')
        temp = self.faces['F'][2]
        self.faces['F'][2] = self.faces['L'][2]
        self.faces['L'][2] = self.faces['B'][2]
        self.faces['B'][2] = self.faces['R'][2]
        self.faces['R'][2] = temp

def solve():
    # Initial state of the cube based on the problem description
    # F=White, U=Orange, R=Blue, L=Green, B=Yellow, D=Red
    initial_state = {
        'F': [['R', 'Y', 'Y'], ['G', 'W', 'W'], ['G', 'O', 'O']], # White
        'U': [['R', 'Y', 'W'], ['B', 'O', 'Y'], ['Y', 'R', 'O']], # Orange
        'R': [['G', 'G', 'G'], ['R', 'B', 'B'], ['B', 'B', 'O']], # Blue
        'L': [['R', 'W', 'W'], ['R', 'G', 'O'], ['W', 'W', 'B']], # Green
        'B': [['Y', 'W', 'R'], ['B', 'Y', 'O'], ['G', 'O', 'B']], # Yellow
        'D': [['B', 'Y', 'Y'], ['R', 'R', 'G'], ['W', 'G', 'O']]  # Red
    }

    # Create a cube instance
    cube = RubiksCube(initial_state)

    # Execute the 5-step algorithm
    cube.R()
    cube.U()
    cube.F()
    cube.L_prime()
    cube.D()

    # Get the final state of the white face (Front)
    final_white_face = cube.faces['F']
    
    # Print the final matrix for the white face
    print(final_white_face)

solve()
<<<E>>>