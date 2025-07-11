import copy

class Cube:
    """Represents a 3x3 Rubik's Cube and its state."""

    def __init__(self, initial_state):
        """Initializes the cube with a given state."""
        # The state is a dictionary mapping face names to 3x3 matrices.
        # F=Front, U=Up, R=Right, B=Back, L=Left, D=Down
        self.faces = copy.deepcopy(initial_state)

    def print_final_face(self, face_key='F'):
        """Prints the final state of the specified face."""
        final_face = self.faces[face_key]
        # The prompt requires printing each element of the final matrix.
        print("The final state of the white face is:")
        print(f"[[{final_face[0][0]},{final_face[0][1]},{final_face[0][2]}],[{final_face[1][0]},{final_face[1][1]},{final_face[1][2]}],[{final_face[2][0]},{final_face[2][1]},{final_face[2][2]}]]")


    def rotate_face_cw(self, face_key):
        """Rotates a face 90 degrees clockwise."""
        face = self.faces[face_key]
        self.faces[face_key] = [
            [face[2][0], face[1][0], face[0][0]],
            [face[2][1], face[1][1], face[0][1]],
            [face[2][2], face[1][2], face[0][2]]
        ]

    def rotate_face_ccw(self, face_key):
        """Rotates a face 90 degrees counter-clockwise."""
        face = self.faces[face_key]
        self.faces[face_key] = [
            [face[0][2], face[1][2], face[2][2]],
            [face[0][1], face[1][1], face[2][1]],
            [face[0][0], face[1][0], face[2][0]]
        ]

    # --- Move Implementations ---
    # Each move rotates the primary face and swaps the adjacent side strips.
    
    def move_R(self): # Clockwise Right
        self.rotate_face_cw('R')
        temp_faces = copy.deepcopy(self.faces)
        for i in range(3):
            self.faces['F'][i][2] = temp_faces['D'][i][2]
            self.faces['U'][i][2] = temp_faces['F'][i][2]
            self.faces['D'][i][2] = temp_faces['B'][2-i][0]
            self.faces['B'][2-i][0] = temp_faces['U'][i][2]

    def move_U(self): # Clockwise Up
        self.rotate_face_cw('U')
        temp_faces = copy.deepcopy(self.faces)
        self.faces['F'][0] = temp_faces['R'][0]
        self.faces['R'][0] = temp_faces['B'][0]
        self.faces['B'][0] = temp_faces['L'][0]
        self.faces['L'][0] = temp_faces['F'][0]

    def move_F(self): # Clockwise Front
        self.rotate_face_cw('F')
        temp_faces = copy.deepcopy(self.faces)
        # Piece-by-piece swaps for clarity
        # Corners: U20->R00->D02->L22->U20 and U22->R20->D00->L02->U22
        # Edges: U21->R10->D01->L12->U21
        self.faces['U'][2][0] = temp_faces['L'][2][2]
        self.faces['U'][2][1] = temp_faces['L'][1][2]
        self.faces['U'][2][2] = temp_faces['L'][0][2]

        self.faces['L'][0][2] = temp_faces['D'][0][0]
        self.faces['L'][1][2] = temp_faces['D'][0][1]
        self.faces['L'][2][2] = temp_faces['D'][0][2]
        
        self.faces['D'][0][0] = temp_faces['R'][2][0]
        self.faces['D'][0][1] = temp_faces['R'][1][0]
        self.faces['D'][0][2] = temp_faces['R'][0][0]

        self.faces['R'][0][0] = temp_faces['U'][2][0]
        self.faces['R'][1][0] = temp_faces['U'][2][1]
        self.faces['R'][2][0] = temp_faces['U'][2][2]

    def move_L_ccw(self): # Counter-Clockwise Left
        self.rotate_face_ccw('L')
        temp_faces = copy.deepcopy(self.faces)
        for i in range(3):
            self.faces['F'][i][0] = temp_faces['D'][i][0]
            self.faces['D'][i][0] = temp_faces['B'][2-i][2]
            self.faces['B'][2-i][2] = temp_faces['U'][i][0]
            self.faces['U'][i][0] = temp_faces['F'][i][0]

    def move_D(self): # Clockwise Down
        self.rotate_face_cw('D')
        temp_faces = copy.deepcopy(self.faces)
        self.faces['F'][2] = temp_faces['L'][2]
        self.faces['L'][2] = temp_faces['B'][2]
        self.faces['B'][2] = temp_faces['R'][2]
        self.faces['R'][2] = temp_faces['F'][2]


def solve():
    """Sets up the cube, runs the algorithm, and prints the result."""
    
    # Map colors to standard faces based on the problem's starting orientation
    # Front (F) is White, Up (U) is Orange, Right (R) is Blue, etc.
    initial_state = {
        'F': [['R','Y','Y'],['G','W','W'],['G','O','O']], # White face
        'U': [['R','Y','W'],['B','O','Y'],['Y','R','O']], # Orange face
        'R': [['G','G','G'],['R','B','B'],['B','B','O']], # Blue face
        'B': [['Y','W','R'],['B','Y','O'],['G','O','B']], # Yellow face
        'L': [['R','W','W'],['R','G','O'],['W','W','B']], # Green face
        'D': [['B','Y','Y'],['R','R','G'],['W','G','O']]  # Red face
    }

    # Create a cube instance
    cube = Cube(initial_state)

    # Execute the 5-step algorithm
    # 1. R
    cube.move_R()
    # 2. U
    cube.move_U()
    # 3. F
    cube.move_F()
    # 4. L'
    cube.move_L_ccw()
    # 5. D
    cube.move_D()

    # Print the final state of the Front (White) face
    cube.print_final_face('F')

solve()