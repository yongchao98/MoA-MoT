import copy

class Cube:
    """
    Represents a Rubik's cube and its transformations.
    The cube is oriented with White as Front (F), Orange as Up (U),
    Blue as Right (R), Green as Left (L), Yellow as Back (B), and Red as Down (D).
    """

    def __init__(self, state):
        """Initializes the cube with a given state."""
        self.faces = copy.deepcopy(state)

    def _rotate_face_cw(self, face_key):
        """Rotates a single face 90 degrees clockwise."""
        face = self.faces[face_key]
        self.faces[face_key] = [list(row) for row in zip(*face[::-1])]

    def _rotate_face_ccw(self, face_key):
        """Rotates a single face 90 degrees counter-clockwise."""
        face = self.faces[face_key]
        self.faces[face_key] = [list(row) for row in zip(*face)][::-1]

    def R(self):
        """Performs a clockwise R (Right) move."""
        self._rotate_face_cw('R')
        # The adjacent faces cycle F -> U -> B -> D -> F
        temp = [row[2] for row in self.faces['F']]
        # D -> F
        for i in range(3): self.faces['F'][i][2] = self.faces['D'][i][2]
        # B(reversed) -> D
        for i in range(3): self.faces['D'][i][2] = self.faces['B'][2-i][0]
        # U -> B(reversed)
        for i in range(3): self.faces['B'][2-i][0] = self.faces['U'][i][2]
        # temp(F) -> U
        for i in range(3): self.faces['U'][i][2] = temp[i]

    def U(self):
        """Performs a clockwise U (Up) move."""
        self._rotate_face_cw('U')
        # The adjacent faces cycle F -> R -> B -> L -> F
        temp = list(self.faces['F'][0])
        self.faces['F'][0] = self.faces['R'][0]
        self.faces['R'][0] = self.faces['B'][0]
        self.faces['B'][0] = self.faces['L'][0]
        self.faces['L'][0] = temp

    def F(self):
        """Performs a clockwise F (Front) move."""
        self._rotate_face_cw('F')
        # The adjacent faces cycle U -> R -> D -> L -> U
        temp_u_bottom = list(self.faces['U'][2])
        # L right col (rev) -> U bottom row
        for i in range(3): self.faces['U'][2][i] = self.faces['L'][2-i][2]
        # D top row (rev) -> L right col
        temp_d_top = list(self.faces['D'][0])
        for i in range(3): self.faces['L'][i][2] = temp_d_top[2-i]
        # R left col (rev) -> D top row
        temp_r_left = [row[0] for row in self.faces['R']]
        for i in range(3): self.faces['D'][0][i] = temp_r_left[2-i]
        # Original U bottom row -> R left col
        for i in range(3): self.faces['R'][i][0] = temp_u_bottom[i]

    def L_prime(self):
        """Performs a counter-clockwise L' (Left) move."""
        self._rotate_face_ccw('L')
        # The adjacent faces cycle F -> U -> B -> D -> F
        temp = [row[0] for row in self.faces['F']]
        self.faces['F'][0] = self.faces['U'][0]
        for i in range(3): self.faces['U'][i][0] = self.faces['B'][2-i][2]
        for i in range(3): self.faces['B'][2-i][2] = self.faces['D'][i][0]
        for i in range(3): self.faces['D'][i][0] = temp[i]
        
    def D(self):
        """Performs a clockwise D (Down) move."""
        self._rotate_face_cw('D')
        # The adjacent faces cycle F -> L -> B -> R -> F
        temp = list(self.faces['F'][2])
        self.faces['F'][2] = self.faces['L'][2]
        self.faces['L'][2] = self.faces['B'][2]
        self.faces['B'][2] = self.faces['R'][2]
        self.faces['R'][2] = temp

def solve_cube_problem():
    """Sets up the cube, runs the algorithm, and prints the result."""
    # Define the jumbled state of the cube based on the problem description
    initial_state = {
        'F': [['R', 'Y', 'Y'], ['G', 'W', 'W'], ['G', 'O', 'O']],  # White face
        'U': [['R', 'Y', 'W'], ['B', 'O', 'Y'], ['Y', 'R', 'O']],  # Orange face
        'R': [['G', 'G', 'G'], ['R', 'B', 'B'], ['B', 'B', 'O']],  # Blue face
        'B': [['Y', 'W', 'R'], ['B', 'Y', 'O'], ['G', 'O', 'B']],  # Yellow face
        'L': [['R', 'W', 'W'], ['R', 'G', 'O'], ['W', 'W', 'B']],  # Green face
        'D': [['B', 'Y', 'Y'], ['R', 'R', 'G'], ['W', 'G', 'O']]   # Red face
    }

    # Create a Cube instance
    cube = Cube(initial_state)

    # Apply the 5-step algorithm: R, U, F, L', D
    cube.R()
    cube.U()
    cube.F()
    cube.L_prime()
    cube.D()

    # The problem asks for the final state of the white face (F)
    final_face = cube.faces['F']
    print(final_face)

solve_cube_problem()
<<<E>>>