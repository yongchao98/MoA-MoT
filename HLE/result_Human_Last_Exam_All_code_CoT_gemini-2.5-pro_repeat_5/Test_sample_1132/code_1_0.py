import copy

class RubiksCube:
    """
    A class to represent and manipulate a Rubik's Cube.
    The cube orientation is fixed: F(ront), U(p), R(ight), L(eft), B(ack), D(own).
    """

    def __init__(self):
        """Initializes the cube to the jumbled state specified in the problem."""
        # Mapping faces to standard orientation:
        # F=Front(White), U=Up(Orange), R=Right(Blue), L=Left(Green), B=Back(Yellow), D=Down(Red)
        self.faces = {
            'F': [['R', 'Y', 'Y'], ['G', 'W', 'W'], ['G', 'O', 'O']],
            'U': [['R', 'Y', 'W'], ['B', 'O', 'Y'], ['Y', 'R', 'O']],
            'R': [['G', 'G', 'G'], ['R', 'B', 'B'], ['B', 'B', 'O']],
            'B': [['Y', 'W', 'R'], ['B', 'Y', 'O'], ['G', 'O', 'B']],
            'L': [['R', 'W', 'W'], ['R', 'G', 'O'], ['W', 'W', 'B']],
            'D': [['B', 'Y', 'Y'], ['R', 'R', 'G'], ['W', 'G', 'O']],
        }

    def _rotate_face_cw(self, face):
        """Rotates a 3x3 face matrix 90 degrees clockwise."""
        return [list(row) for row in zip(*face[::-1])]

    def _rotate_face_ccw(self, face):
        """Rotates a 3x3 face matrix 90 degrees counter-clockwise."""
        return [list(row) for row in zip(*face)][::-1]

    def move_R(self):
        """Performs a clockwise R (Right) move."""
        old_faces = copy.deepcopy(self.faces)
        self.faces['R'] = self._rotate_face_cw(old_faces['R'])
        # Adjacency cycle: F -> U -> B -> D -> F (with twists for B)
        for i in range(3):
            self.faces['F'][i][2] = old_faces['D'][i][2]
            self.faces['U'][i][2] = old_faces['F'][i][2]
            self.faces['B'][i][0] = old_faces['U'][2 - i][2]
            self.faces['D'][i][2] = old_faces['B'][2 - i][0]

    def move_U(self):
        """Performs a clockwise U (Up) move."""
        old_faces = copy.deepcopy(self.faces)
        self.faces['U'] = self._rotate_face_cw(old_faces['U'])
        # Adjacency cycle: F -> R -> B -> L -> F
        self.faces['F'][0] = old_faces['R'][0]
        self.faces['R'][0] = old_faces['B'][0]
        self.faces['B'][0] = old_faces['L'][0]
        self.faces['L'][0] = old_faces['F'][0]

    def move_F(self):
        """Performs a clockwise F (Front) move."""
        old_faces = copy.deepcopy(self.faces)
        self.faces['F'] = self._rotate_face_cw(old_faces['F'])
        # Adjacency cycle: U -> R -> D -> L -> U
        for i in range(3):
            self.faces['U'][2][i] = old_faces['L'][i][2]
            self.faces['R'][i][0] = old_faces['U'][2][i]
            self.faces['D'][0][i] = old_faces['R'][i][0]
            self.faces['L'][i][2] = old_faces['D'][0][i]

    def move_L_prime(self):
        """Performs a counter-clockwise L' (Left) move."""
        old_faces = copy.deepcopy(self.faces)
        self.faces['L'] = self._rotate_face_ccw(old_faces['L'])
        # Adjacency cycle (L'): F -> D -> B -> U -> F (with twists for B)
        for i in range(3):
            self.faces['F'][i][0] = old_faces['U'][i][0]
            self.faces['U'][i][0] = old_faces['B'][2 - i][2]
            self.faces['B'][2 - i][2] = old_faces['D'][i][0]
            self.faces['D'][i][0] = old_faces['F'][i][0]

    def move_D(self):
        """Performs a clockwise D (Down) move."""
        old_faces = copy.deepcopy(self.faces)
        self.faces['D'] = self._rotate_face_cw(old_faces['D'])
        # Adjacency cycle: F -> L -> B -> R -> F
        self.faces['F'][2] = old_faces['L'][2]
        self.faces['L'][2] = old_faces['B'][2]
        self.faces['B'][2] = old_faces['R'][2]
        self.faces['R'][2] = old_faces['F'][2]

    def get_face(self, face_name):
        """Returns the state of a specific face."""
        return self.faces[face_name]

# Main execution block
if __name__ == "__main__":
    # Create a cube instance with the specified jumbled state
    cube = RubiksCube()

    # Apply the 5-step algorithm
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

    # Get the final state of the white face (Front face)
    final_white_face = cube.get_face('F')

    # Print the result to compare with answer choices
    print(final_white_face)
<<<F>>>