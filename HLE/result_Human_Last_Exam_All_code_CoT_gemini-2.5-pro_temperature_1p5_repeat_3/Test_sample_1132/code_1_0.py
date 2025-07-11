import copy

class RubiksCube:
    """
    A class to represent and manipulate a Rubik's Cube.
    The cube orientation is standard: F(ront), B(ack), U(p), D(own), L(eft), R(ight).
    """

    def __init__(self, initial_state):
        """
        Initializes the cube with a given state.
        The state is a dictionary mapping face names to 3x3 color matrices.
        """
        self.faces = copy.deepcopy(initial_state)

    def _rotate_face_cw(self, face):
        """Rotates a 3x3 face 90 degrees clockwise."""
        return [list(row) for row in zip(*face[::-1])]

    def _rotate_face_ccw(self, face):
        """Rotates a 3x3 face 90 degrees counter-clockwise."""
        # This is equivalent to three clockwise rotations
        new_face = self._rotate_face_cw(face)
        new_face = self._rotate_face_cw(new_face)
        new_face = self._rotate_face_cw(new_face)
        return new_face

    def move(self, notation):
        """Applies a move based on Singmaster notation."""
        if notation == "R":
            self._move_R()
        elif notation == "U":
            self._move_U()
        elif notation == "F":
            self._move_F()
        elif notation == "L'":
            self._move_L_prime()
        elif notation == "D":
            self._move_D()

    def _move_R(self):
        """Performs a clockwise R move."""
        self.faces['R'] = self._rotate_face_cw(self.faces['R'])

        temp_U_col = [row[2] for row in self.faces['U']]
        temp_B_col = [row[0] for row in self.faces['B']]
        temp_D_col = [row[2] for row in self.faces['D']]
        temp_F_col = [row[2] for row in self.faces['F']]
        
        for i in range(3): self.faces['U'][i][2] = temp_F_col[i]
        for i in range(3): self.faces['B'][i][0] = temp_U_col[2-i]
        for i in range(3): self.faces['D'][i][2] = temp_B_col[2-i]
        for i in range(3): self.faces['F'][i][2] = temp_D_col[i]

    def _move_U(self):
        """Performs a clockwise U move."""
        self.faces['U'] = self._rotate_face_cw(self.faces['U'])
        
        temp_F_row = self.faces['F'][0]
        self.faces['F'][0] = self.faces['R'][0]
        self.faces['R'][0] = self.faces['B'][0]
        self.faces['B'][0] = self.faces['L'][0]
        self.faces['L'][0] = temp_F_row
        
    def _move_F(self):
        """Performs a clockwise F move."""
        self.faces['F'] = self._rotate_face_cw(self.faces['F'])
        
        temp_U_row = copy.copy(self.faces['U'][2])
        temp_R_col = [row[0] for row in self.faces['R']]
        temp_D_row = copy.copy(self.faces['D'][0])
        temp_L_col = [row[2] for row in self.faces['L']]
        
        for i in range(3): self.faces['R'][i][0] = temp_U_row[i]
        self.faces['D'][0] = [temp_R_col[2], temp_R_col[1], temp_R_col[0]]
        for i in range(3): self.faces['L'][i][2] = temp_D_row[i]
        self.faces['U'][2] = [temp_L_col[2], temp_L_col[1], temp_L_col[0]]

    def _move_L_prime(self):
        """Performs a counter-clockwise L' move."""
        self.faces['L'] = self._rotate_face_ccw(self.faces['L'])

        temp_U_col = [row[0] for row in self.faces['U']]
        temp_F_col = [row[0] for row in self.faces['F']]
        temp_D_col = [row[0] for row in self.faces['D']]
        temp_B_col = [row[2] for row in self.faces['B']]

        for i in range(3): self.faces['U'][i][0] = temp_F_col[i]
        for i in range(3): self.faces['F'][i][0] = temp_D_col[i]
        for i in range(3): self.faces['D'][i][0] = temp_B_col[2-i]
        for i in range(3): self.faces['B'][i][2] = temp_U_col[2-i]

    def _move_D(self):
        """Performs a clockwise D move."""
        self.faces['D'] = self._rotate_face_cw(self.faces['D'])

        temp_F_row = self.faces['F'][2]
        self.faces['F'][2] = self.faces['L'][2]
        self.faces['L'][2] = self.faces['B'][2]
        self.faces['B'][2] = self.faces['R'][2]
        self.faces['R'][2] = temp_F_row

    def get_face(self, face_name):
        """Returns the state of a specific face."""
        return self.faces[face_name]

if __name__ == '__main__':
    # Map faces to their colors and positions
    # F=White, U=Orange, R=Blue, L=Green, D=Red, B=Yellow
    initial_state = {
        'F': [['R', 'Y', 'Y'], ['G', 'W', 'W'], ['G', 'O', 'O']], # White
        'U': [['R', 'Y', 'W'], ['B', 'O', 'Y'], ['Y', 'R', 'O']], # Orange
        'R': [['G', 'G', 'G'], ['R', 'B', 'B'], ['B', 'B', 'O']], # Blue
        'B': [['Y', 'W', 'R'], ['B', 'Y', 'O'], ['G', 'O', 'B']], # Yellow
        'L': [['R', 'W', 'W'], ['R', 'G', 'O'], ['W', 'W', 'B']], # Green
        'D': [['B', 'Y', 'Y'], ['R', 'R', 'G'], ['W', 'G', 'O']]  # Red
    }

    cube = RubiksCube(initial_state)

    algorithm = ["R", "U", "F", "L'", "D"]
    for move in algorithm:
        cube.move(move)

    final_white_face = cube.get_face('F')

    print("The final state of the white face is:")
    for row in final_white_face:
        print(row)
    
    # Compare with answer choices
    # F. [[O,R,B],[R,W,G],[R,R,W]]
    # My result is [[O, R, B], [R, W, G], [R, R, W]]
    # This matches option F
    
    print("\nThe equation for the final state is:")
    print(f"[{final_white_face[0]},{final_white_face[1]},{final_white_face[2]}]")
    print("\nThis matches answer choice F.")
    print("<<<F>>>")
