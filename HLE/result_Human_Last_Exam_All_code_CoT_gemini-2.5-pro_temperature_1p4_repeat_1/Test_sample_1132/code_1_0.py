import copy

class RubiksCube:
    """A class to represent and manipulate a 3x3 Rubik's Cube."""

    def __init__(self, faces):
        """
        Initializes the cube.
        Faces is a dictionary {'U', 'D', 'L', 'R', 'F', 'B'} where each
        value is a 3x3 matrix (list of lists) representing the colors.
        """
        self.faces = faces

    def get_face(self, face_name):
        return self.faces[face_name]

    def _rotate_face_clockwise(self, face_name):
        """Rotates a given face 90 degrees clockwise."""
        face = self.faces[face_name]
        self.faces[face_name] = [
            [face[2][0], face[1][0], face[0][0]],
            [face[2][1], face[1][1], face[0][1]],
            [face[2][2], face[1][2], face[0][2]],
        ]

    def _rotate_face_counter_clockwise(self, face_name):
        """Rotates a given face 90 degrees counter-clockwise."""
        face = self.faces[face_name]
        self.faces[face_name] = [
            [face[0][2], face[1][2], face[2][2]],
            [face[0][1], face[1][1], face[2][1]],
            [face[0][0], face[1][0], face[2][0]],
        ]

    def R(self):
        """Performs a clockwise R (Right) move."""
        old_faces = copy.deepcopy(self.faces)
        self._rotate_face_clockwise('R')
        # Cycle strips: F(right col) -> U(right col) -> B(left col, reversed) -> D(right col) -> F(right col)
        for i in range(3):
            self.faces['U'][i][2] = old_faces['F'][i][2]
            self.faces['B'][i][0] = old_faces['U'][2 - i][2]
            self.faces['D'][i][2] = old_faces['B'][2 - i][0]
            self.faces['F'][i][2] = old_faces['D'][i][2]
        return self

    def U(self):
        """Performs a clockwise U (Up) move."""
        old_faces = copy.deepcopy(self.faces)
        self._rotate_face_clockwise('U')
        # Cycle strips: F(top row) -> R(top row) -> B(top row) -> L(top row) -> F(top row)
        self.faces['F'][0] = old_faces['R'][0]
        self.faces['R'][0] = old_faces['B'][0]
        self.faces['B'][0] = old_faces['L'][0]
        self.faces['L'][0] = old_faces['F'][0]
        return self

    def F(self):
        """Performs a clockwise F (Front) move."""
        old_faces = copy.deepcopy(self.faces)
        self._rotate_face_clockwise('F')
        # Cycle strips: U(bottom row) -> R(left col) -> D(top row, reversed) -> L(right col) -> U(bottom row, reversed)
        for i in range(3):
            self.faces['R'][i][0] = old_faces['U'][2][i]
            self.faces['D'][0][i] = old_faces['R'][2 - i][0]
            self.faces['L'][i][2] = old_faces['D'][0][i]
            self.faces['U'][2][i] = old_faces['L'][2 - i][2]
        return self

    def L_prime(self):
        """Performs a counter-clockwise L' (Left) move."""
        old_faces = copy.deepcopy(self.faces)
        self._rotate_face_counter_clockwise('L')
        # Cycle strips: F(left col) -> D(left col) -> B(right col, reversed) -> U(left col) -> F(left col)
        for i in range(3):
            self.faces['D'][i][0] = old_faces['F'][i][0]
            self.faces['B'][i][2] = old_faces['D'][2 - i][0]
            self.faces['U'][i][0] = old_faces['B'][2 - i][2]
            self.faces['F'][i][0] = old_faces['U'][i][0]
        return self

    def D(self):
        """Performs a clockwise D (Down) move."""
        old_faces = copy.deepcopy(self.faces)
        self._rotate_face_clockwise('D')
        # Cycle strips: F(bottom row) -> L(bottom row) -> B(bottom row) -> R(bottom row) -> F(bottom row)
        self.faces['F'][2] = old_faces['L'][2]
        self.faces['L'][2] = old_faces['B'][2]
        self.faces['B'][2] = old_faces['R'][2]
        self.faces['R'][2] = old_faces['F'][2]
        return self

if __name__ == '__main__':
    # Define the initial jumbled state of the cube
    # F=White, U=Orange, R=Blue, B=Yellow, L=Green, D=Red
    initial_faces = {
        'F': [['R', 'Y', 'Y'], ['G', 'W', 'W'], ['G', 'O', 'O']],  # White face
        'U': [['R', 'Y', 'W'], ['B', 'O', 'Y'], ['Y', 'R', 'O']],  # Orange face
        'R': [['G', 'G', 'G'], ['R', 'B', 'B'], ['B', 'B', 'O']],  # Blue face
        'B': [['Y', 'W', 'R'], ['B', 'Y', 'O'], ['G', 'O', 'B']],  # Yellow face
        'L': [['R', 'W', 'W'], ['R', 'G', 'O'], ['W', 'W', 'B']],  # Green face
        'D': [['B', 'Y', 'Y'], ['R', 'R', 'G'], ['W', 'G', 'O']]   # Red face
    }

    # Create a cube instance
    cube = RubiksCube(initial_faces)

    # Perform the 5-step algorithm
    cube.R()
    cube.U()
    cube.F()
    cube.L_prime()
    cube.D()

    # Get the final state of the white (Front) face
    final_white_face = cube.get_face('F')

    # Print the result
    print("The final state of the white face is:")
    for row in final_white_face:
        print(row)
<<<E>>>