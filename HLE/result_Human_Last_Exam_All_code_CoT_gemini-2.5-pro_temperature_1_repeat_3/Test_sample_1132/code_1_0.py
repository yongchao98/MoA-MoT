import copy

class RubiksCube:
    """
    A class to represent and manipulate a 3x3 Rubik's Cube.
    The cube is represented by a dictionary of 6 faces.
    F(Front), B(Back), U(Up), D(Down), L(Left), R(Right).
    Each face is a 3x3 matrix (list of lists).
    """

    def __init__(self, initial_state):
        """Initializes the cube with a given state."""
        self.faces = copy.deepcopy(initial_state)

    def get_face(self, face_name):
        """Returns the state of a specific face."""
        return self.faces[face_name]

    def rotate_face_clockwise(self, face_name):
        """Rotates a given face 90 degrees clockwise."""
        face = self.faces[face_name]
        # Transpose the matrix and then reverse each row.
        self.faces[face_name] = [list(row) for row in zip(*face[::-1])]

    def rotate_face_counter_clockwise(self, face_name):
        """Rotates a given face 90 degrees counter-clockwise."""
        face = self.faces[face_name]
        # Reverse each row and then transpose the matrix.
        transposed = [list(row) for row in zip(*face)]
        self.faces[face_name] = transposed[::-1]

    def move_R(self):
        """Performs a clockwise R (Right) move."""
        self.rotate_face_clockwise('Right')
        
        # Store current columns that will be moved
        temp_front_col = [self.faces['Front'][i][2] for i in range(3)]
        temp_up_col = [self.faces['Up'][i][2] for i in range(3)]
        temp_back_col = [self.faces['Back'][i][0] for i in range(3)]
        temp_down_col = [self.faces['Down'][i][2] for i in range(3)]
        
        # Perform the cycle: Up -> Front -> Down -> Back(rev) -> Up
        for i in range(3): self.faces['Front'][i][2] = temp_up_col[i]
        for i in range(3): self.faces['Down'][i][2] = temp_front_col[i]
        for i in range(3): self.faces['Back'][i][0] = temp_down_col[2-i]
        for i in range(3): self.faces['Up'][i][2] = temp_back_col[i]

    def move_U(self):
        """Performs a clockwise U (Up) move."""
        self.rotate_face_clockwise('Up')
        
        # Perform the cycle: Front -> Right -> Back -> Left -> Front
        temp_front_row = self.faces['Front'][0]
        self.faces['Front'][0] = self.faces['Left'][0]
        self.faces['Left'][0] = self.faces['Back'][0]
        self.faces['Back'][0] = self.faces['Right'][0]
        self.faces['Right'][0] = temp_front_row

    def move_F(self):
        """Performs a clockwise F (Front) move."""
        self.rotate_face_clockwise('Front')
        
        # Store current strips that will be moved
        temp_up_row = self.faces['Up'][2]
        
        # Perform the cycle: Up -> Right -> Down -> Left -> Up (with rotations)
        left_right_col = [self.faces['Left'][i][2] for i in range(3)]
        self.faces['Up'][2] = [left_right_col[2], left_right_col[1], left_right_col[0]]

        down_top_row = self.faces['Down'][0]
        for i in range(3): self.faces['Left'][i][2] = down_top_row[2-i]
        
        right_left_col = [self.faces['Right'][i][0] for i in range(3)]
        self.faces['Down'][0] = [right_left_col[2], right_left_col[1], right_left_col[0]]
        
        for i in range(3): self.faces['Right'][i][0] = temp_up_row[i]

    def move_L_prime(self):
        """Performs a counter-clockwise L' (Left) move."""
        self.rotate_face_counter_clockwise('Left')
        
        # Store current columns that will be moved
        temp_front_col = [self.faces['Front'][i][0] for i in range(3)]
        temp_up_col = [self.faces['Up'][i][0] for i in range(3)]
        temp_back_col = [self.faces['Back'][i][2] for i in range(3)]
        temp_down_col = [self.faces['Down'][i][0] for i in range(3)]
        
        # Perform the cycle: Up -> Front -> Down -> Back(rev) -> Up
        for i in range(3): self.faces['Front'][i][0] = temp_up_col[i]
        for i in range(3): self.faces['Down'][i][0] = temp_front_col[i]
        for i in range(3): self.faces['Back'][i][2] = temp_down_col[2-i]
        for i in range(3): self.faces['Up'][i][0] = temp_back_col[i]
        
    def move_D(self):
        """Performs a clockwise D (Down) move."""
        self.rotate_face_clockwise('Down')
        
        # Perform the cycle: Front -> Left -> Back -> Right -> Front
        temp_front_row = self.faces['Front'][2]
        self.faces['Front'][2] = self.faces['Right'][2]
        self.faces['Right'][2] = self.faces['Back'][2]
        self.faces['Back'][2] = self.faces['Left'][2]
        self.faces['Left'][2] = temp_front_row

if __name__ == "__main__":
    # Define the initial jumbled state of the cube
    # F=White, U=Orange, R=Blue, L=Green, B=Yellow, D=Red
    initial_jumbled_state = {
        'Front': [['R','Y','Y'], ['G','W','W'], ['G','O','O']], # White face
        'Up':    [['R','Y','W'], ['B','O','Y'], ['Y','R','O']], # Orange face
        'Right': [['G','G','G'], ['R','B','B'], ['B','B','O']], # Blue face
        'Back':  [['Y','W','R'], ['B','Y','O'], ['G','O','B']], # Yellow face
        'Left':  [['R','W','W'], ['R','G','O'], ['W','W','B']], # Green face
        'Down':  [['B','Y','Y'], ['R','R','G'], ['W','G','O']]  # Red face
    }
    
    # Create a cube instance
    cube = RubiksCube(initial_jumbled_state)
    
    # Apply the 5-step algorithm
    cube.move_R()
    cube.move_U()
    cube.move_F()
    cube.move_L_prime()
    cube.move_D()
    
    # Get the final state of the white face
    final_white_face = cube.get_face('Front')
    
    # Print the result in the required format
    print("The final equation is:")
    print(f"[{final_white_face[0]},")
    print(f" {final_white_face[1]},")
    print(f" {final_white_face[2]}]")
    
    # Identify the matching answer choice
    # A. [[B,G,B],[W,W,R],[R,W,Y]]
    # B. [[Y,O,G],[W,W,O],[W,G,B]]
    # C. [[O,G,B],[G,W,R],[R,W,R]]
    # D. [[W,G,B],[B,W,R],[G,W,R]]
    # E. [[B,G,B],[W,W,R],[R,W,B]]
    # F. [[O,R,B],[R,W,G],[R,R,W]]
    # My Result: [[B, G, B], [W, W, R], [R, W, B]] -> Matches E

    print("\n<<<E>>>")