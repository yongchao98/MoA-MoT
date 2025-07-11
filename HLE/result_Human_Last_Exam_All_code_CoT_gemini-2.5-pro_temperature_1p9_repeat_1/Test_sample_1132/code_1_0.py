import copy

class Cube:
    """
    A class to represent a Rubik's Cube and its movements.
    The cube state is stored as a dictionary of faces, each being a 3x3 matrix.
    Standard orientation is assumed:
    - front: White face
    - up:    Orange face
    - right: Blue face
    - left:  Green face
    - down:  Red face
    - back:  Yellow face
    """
    def __init__(self, state):
        self.faces = state

    def get_white_face(self):
        return self.faces['front']

    def rotate_face_cw(self, face_name):
        # Performs a 90-degree clockwise rotation on a face's matrix
        face = self.faces[face_name]
        self.faces[face_name] = [[face[2][0], face[1][0], face[0][0]],
                                 [face[2][1], face[1][1], face[0][1]],
                                 [face[2][2], face[1][2], face[0][2]]]

    def rotate_face_ccw(self, face_name):
        # Performs a 90-degree counter-clockwise rotation on a face's matrix
        face = self.faces[face_name]
        self.faces[face_name] = [[face[0][2], face[1][2], face[2][2]],
                                 [face[0][1], face[1][1], face[2][1]],
                                 [face[0][0], face[1][0], face[2][0]]]

    def move(self, notation):
        # A deep copy reads from an unmodified state while writing to the new state.
        old_faces = copy.deepcopy(self.faces)

        if notation == 'R': # Right face clockwise
            self.rotate_face_cw('right')
            # Cycle: Front -> Up -> Back -> Down -> Front
            for i in range(3):
                self.faces['up'][i][2] = old_faces['front'][i][2]
                self.faces['back'][2 - i][0] = old_faces['up'][i][2]
                self.faces['down'][i][2] = old_faces['back'][2 - i][0]
                self.faces['front'][i][2] = old_faces['down'][i][2]

        elif notation == 'U': # Up face clockwise
            self.rotate_face_cw('up')
            # Cycle: Front -> Right -> Back -> Left -> Front
            self.faces['front'][0] = old_faces['right'][0]
            self.faces['right'][0] = old_faces['back'][0]
            self.faces['back'][0] = old_faces['left'][0]
            self.faces['left'][0] = old_faces['front'][0]

        elif notation == 'F': # Front face clockwise
            self.rotate_face_cw('front')
            # Cycle: Up -> Right -> Down -> Left -> Up
            for i in range(3):
                self.faces['right'][i][0] = old_faces['up'][2][i]
                self.faces['down'][0][i] = old_faces['right'][2-i][0]
                self.faces['left'][i][2] = old_faces['down'][0][2-i]
                self.faces['up'][2][i] = old_faces['left'][2-i][2]

        elif notation == "L'": # Left face counter-clockwise
            self.rotate_face_ccw('left')
            # Cycle: Front -> Down -> Back -> Up -> Front
            for i in range(3):
                self.faces['down'][i][0] = old_faces['front'][i][0]
                self.faces['back'][2-i][2] = old_faces['down'][i][0]
                self.faces['up'][i][0] = old_faces['back'][2-i][2]
                self.faces['front'][i][0] = old_faces['up'][i][0]

        elif notation == 'D': # Down face clockwise
            self.rotate_face_cw('down')
            # Cycle: Front -> Left -> Back -> Right -> Front
            self.faces['front'][2] = old_faces['left'][2]
            self.faces['left'][2] = old_faces['back'][2]
            self.faces['back'][2] = old_faces['right'][2]
            self.faces['right'][2] = old_faces['front'][2]

def main():
    # Initial state of the jumbled cube as per the problem description
    initial_state = {
        'front': [['R', 'Y', 'Y'], ['G', 'W', 'W'], ['G', 'O', 'O']],
        'up':    [['R', 'Y', 'W'], ['B', 'O', 'Y'], ['Y', 'R', 'O']],
        'right': [['G', 'G', 'G'], ['R', 'B', 'B'], ['B', 'B', 'O']],
        'back':  [['Y', 'W', 'R'], ['B', 'Y', 'O'], ['G', 'O', 'B']],
        'left':  [['R', 'W', 'W'], ['R', 'G', 'O'], ['W', 'W', 'B']],
        'down':  [['B', 'Y', 'Y'], ['R', 'R', 'G'], ['W', 'G', 'O']]
    }

    # The algorithm to execute
    algorithm = ['R', 'U', 'F', "L'", 'D']

    # Create a cube instance
    my_cube = Cube(initial_state)

    # Apply each move in the algorithm
    for move in algorithm:
        my_cube.move(move)

    # Get the final state of the white face
    final_white_face = my_cube.get_white_face()

    # Format the output to match the answer choices, e.g., [[O,R,B],[R,W,G],[R,R,W]]
    output_str = "[[" + ",".join(final_white_face[0]) + "],[" + ",".join(final_white_face[1]) + "],[" + ",".join(final_white_face[2]) + "]]"
    print(output_str)

if __name__ == "__main__":
    main()