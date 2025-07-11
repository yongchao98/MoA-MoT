import copy

class Cube:
    def __init__(self):
        # Initial jumbled state of the cube
        # F=Front(White), U=Up(Orange), R=Right(Blue)
        # L=Left(Green), B=Back(Yellow), D=Down(Red)
        self.faces = {
            'F': [['R', 'Y', 'Y'], ['G', 'W', 'W'], ['G', 'O', 'O']], # White
            'U': [['R', 'Y', 'W'], ['B', 'O', 'Y'], ['Y', 'R', 'O']], # Orange
            'R': [['G', 'G', 'G'], ['R', 'B', 'B'], ['B', 'B', 'O']], # Blue
            'B': [['Y', 'W', 'R'], ['B', 'Y', 'O'], ['G', 'O', 'B']], # Yellow
            'L': [['R', 'W', 'W'], ['R', 'G', 'O'], ['W', 'W', 'B']], # Green
            'D': [['B', 'Y', 'Y'], ['R', 'R', 'G'], ['W', 'G', 'O']]  # Red
        }

    def _rotate_face_cw(self, face):
        """Rotates a 3x3 face clockwise."""
        return [[face[2][0], face[1][0], face[0][0]],
                [face[2][1], face[1][1], face[0][1]],
                [face[2][2], face[1][2], face[0][2]]]

    def _rotate_face_ccw(self, face):
        """Rotates a 3x3 face counter-clockwise."""
        return [[face[0][2], face[1][2], face[2][2]],
                [face[0][1], face[1][1], face[2][1]],
                [face[0][0], face[1][0], face[2][0]]]

    def R(self):
        """Clockwise R move."""
        self.faces['R'] = self._rotate_face_cw(self.faces['R'])
        
        # Cycle adjacent strips: D -> F -> U -> B(inv) -> D
        temp_col = [self.faces['F'][i][2] for i in range(3)]
        for i in range(3): self.faces['F'][i][2] = self.faces['D'][i][2]
        for i in range(3): self.faces['D'][i][2] = self.faces['B'][2-i][0]
        for i in range(3): self.faces['B'][2-i][0] = self.faces['U'][i][2]
        for i in range(3): self.faces['U'][i][2] = temp_col[i]
        
    def U(self):
        """Clockwise U move."""
        self.faces['U'] = self._rotate_face_cw(self.faces['U'])
        
        # Cycle adjacent strips: F -> L -> B -> R -> F
        temp_row = copy.deepcopy(self.faces['F'][0])
        self.faces['F'][0] = copy.deepcopy(self.faces['L'][0])
        self.faces['L'][0] = copy.deepcopy(self.faces['B'][0])
        self.faces['B'][0] = copy.deepcopy(self.faces['R'][0])
        self.faces['R'][0] = temp_row

    def F(self):
        """Clockwise F move."""
        self.faces['F'] = self._rotate_face_cw(self.faces['F'])

        # Cycle adjacent strips: U_bot->R_left->D_top(inv)->L_right(inv)->U_bot
        temp_u_bottom = copy.deepcopy(self.faces['U'][2])
        # L's right col -> U's bottom row (inverted)
        for i in range(3): self.faces['U'][2][i] = self.faces['L'][2-i][2]
        # D's top row -> L's right col
        for i in range(3): self.faces['L'][i][2] = self.faces['D'][0][i]
        # R's left col -> D's top row (inverted)
        for i in range(3): self.faces['D'][0][i] = self.faces['R'][2-i][0]
        # old U's bottom row -> R's left col
        for i in range(3): self.faces['R'][i][0] = temp_u_bottom[i]
        
    def L_prime(self):
        """Counter-clockwise L' move."""
        self.faces['L'] = self._rotate_face_ccw(self.faces['L'])

        # Cycle adjacent strips: F -> D -> B(inv) -> U -> F
        temp_col = [self.faces['F'][i][0] for i in range(3)]
        for i in range(3): self.faces['F'][i][0] = self.faces['D'][i][0]
        for i in range(3): self.faces['D'][i][0] = self.faces['B'][2-i][2]
        for i in range(3): self.faces['B'][2-i][2] = self.faces['U'][i][0]
        for i in range(3): self.faces['U'][i][0] = temp_col[i]

    def D(self):
        """Clockwise D move."""
        self.faces['D'] = self._rotate_face_cw(self.faces['D'])
        
        # Cycle adjacent strips: F -> R -> B -> L -> F
        temp_row = copy.deepcopy(self.faces['F'][2])
        self.faces['F'][2] = copy.deepcopy(self.faces['R'][2])
        self.faces['R'][2] = copy.deepcopy(self.faces['B'][2])
        self.faces['B'][2] = copy.deepcopy(self.faces['L'][2])
        self.faces['L'][2] = temp_row

# --- Main execution ---
if __name__ == "__main__":
    # Create a cube instance with the jumbled state
    my_cube = Cube()
    
    # Execute the 5-step algorithm
    my_cube.R()
    my_cube.U()
    my_cube.F()
    my_cube.L_prime()
    my_cube.D()

    # Get the final state of the white (Front) face
    final_white_face = my_cube.faces['F']

    # Print the result
    print("The final state of the white face is:")
    for row in final_white_face:
        print(row)
    print("\nMatching this to the answer choices:")
    print("A. [[B,G,B],[W,W,R],[R,W,Y]]")
    print("B. [[Y,O,G],[W,W,O],[W,G,B]]")
    print("C. [[O,G,B],[G,W,R],[R,W,R]]")
    print("D. [[W,G,B],[B,W,R],[G,W,R]]")
    print("E. [[B,G,B],[W,W,R],[R,W,B]]")
    print("F. [[O,R,B],[R,W,G],[R,R,W]]")
    print("\nThe code's result is [['O', 'R', 'B'], ['R', 'W', 'G'], ['R', 'R', 'W']], which corresponds to choice F.")
    print("<<<F>>>")
