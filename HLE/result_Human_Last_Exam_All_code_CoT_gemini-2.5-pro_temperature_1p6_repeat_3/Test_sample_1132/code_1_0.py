import copy

class RubiksCube:
    def __init__(self, initial_state):
        self.faces = copy.deepcopy(initial_state)

    def _rotate_face_cw(self, face_name):
        face = self.faces[face_name]
        self.faces[face_name] = [[face[2 - j][i] for j in range(3)] for i in range(3)]

    def _rotate_face_ccw(self, face_name):
        face = self.faces[face_name]
        self.faces[face_name] = [[face[j][2 - i] for j in range(3)] for i in range(3)]

    def R(self):
        self._rotate_face_cw('R')
        temp_f = [self.faces['F'][i][2] for i in range(3)]
        temp_u = [self.faces['U'][i][2] for i in range(3)]
        temp_b = [self.faces['B'][i][0] for i in range(3)]
        temp_d = [self.faces['D'][i][2] for i in range(3)]
        for i in range(3): self.faces['U'][i][2] = temp_f[i]
        for i in range(3): self.faces['B'][2 - i][0] = temp_u[i]
        for i in range(3): self.faces['D'][i][2] = temp_b[2-i]
        for i in range(3): self.faces['F'][i][2] = temp_d[i]

    def U(self):
        self._rotate_face_cw('U')
        temp_f = self.faces['F'][0][:]
        temp_r = self.faces['R'][0][:]
        temp_b = self.faces['B'][0][:]
        temp_l = self.faces['L'][0][:]
        self.faces['F'][0] = temp_r
        self.faces['R'][0] = temp_b
        self.faces['B'][0] = temp_l
        self.faces['L'][0] = temp_f

    def F(self):
        self._rotate_face_cw('F')
        temp_u = self.faces['U'][2][:]
        temp_r = [self.faces['R'][i][0] for i in range(3)]
        temp_d = self.faces['D'][0][:]
        temp_l = [self.faces['L'][i][2] for i in range(3)]
        for i in range(3): self.faces['R'][i][0] = temp_u[i]
        for i in range(3): self.faces['D'][0][i] = temp_r[2 - i]
        for i in range(3): self.faces['L'][i][2] = temp_d[i]
        for i in range(3): self.faces['U'][2][i] = temp_l[2 - i]

    def L_prime(self):
        self._rotate_face_ccw('L')
        temp_f = [self.faces['F'][i][0] for i in range(3)]
        temp_u = [self.faces['U'][i][0] for i in range(3)]
        temp_b = [self.faces['B'][i][2] for i in range(3)]
        temp_d = [self.faces['D'][i][0] for i in range(3)]
        for i in range(3): self.faces['U'][i][0] = temp_f[i]
        for i in range(3): self.faces['B'][2-i][2] = temp_u[i]
        for i in range(3): self.faces['D'][i][0] = temp_b[2-i]
        for i in range(3): self.faces['F'][i][0] = temp_d[i]
        
    def D(self):
        self._rotate_face_cw('D')
        temp_f = self.faces['F'][2][:]
        temp_r = self.faces['R'][2][:]
        temp_b = self.faces['B'][2][:]
        temp_l = self.faces['L'][2][:]
        self.faces['F'][2] = temp_l
        self.faces['L'][2] = temp_b
        self.faces['B'][2] = temp_r
        self.faces['R'][2] = temp_f

if __name__ == '__main__':
    # Define faces based on problem statement: White=Front, Orange=Up, Blue=Right
    # This implies Green=Left, Red=Bottom, Yellow=Back
    initial_cube_state = {
        'F': [['R', 'Y', 'Y'], ['G', 'W', 'W'], ['G', 'O', 'O']],  # White face
        'U': [['R', 'Y', 'W'], ['B', 'O', 'Y'], ['Y', 'R', 'O']],  # Orange face
        'R': [['G', 'G', 'G'], ['R', 'B', 'B'], ['B', 'B', 'O']],  # Blue face
        'B': [['Y', 'W', 'R'], ['B', 'Y', 'O'], ['G', 'O', 'B']],  # Yellow face
        'L': [['R', 'W', 'W'], ['R', 'G', 'O'], ['W', 'W', 'B']],  # Green face
        'D': [['B', 'Y', 'Y'], ['R', 'R', 'G'], ['W', 'G', 'O']],  # Red face
    }

    cube = RubiksCube(initial_cube_state)
    
    # Apply the 5-step algorithm
    cube.R()
    cube.U()
    cube.F()
    cube.L_prime()
    cube.D()

    final_white_face = cube.faces['F']
    
    # Print the result in the required format
    print("The final state of the white face is:")
    for row in final_white_face:
        print(row)
    
    print("\nThis corresponds to the equation:")
    row1 = " ".join(final_white_face[0])
    row2 = " ".join(final_white_face[1])
    row3 = " ".join(final_white_face[2])
    print(f"[[{row1}], [{row2}], [{row3}]]")

    # This part is just to help the user identify the correct option
    # It will not be part of the thinking process, but compares the result to the options.
    options = {
        "A": [['B', 'G', 'B'], ['W', 'W', 'R'], ['R', 'W', 'Y']],
        "B": [['Y', 'O', 'G'], ['W', 'W', 'O'], ['W', 'G', 'B']],
        "C": [['O', 'G', 'B'], ['G', 'W', 'R'], ['R', 'W', 'R']],
        "D": [['W', 'G', 'B'], ['B', 'W', 'R'], ['G', 'W', 'R']],
        "E": [['B', 'G', 'B'], ['W', 'W', 'R'], ['R', 'W', 'B']],
        "F": [['O', 'R', 'B'], ['R', 'W', 'G'], ['R', 'R', 'W']]
    }

    for letter, matrix in options.items():
        if final_white_face == matrix:
            print(f"\nThis matches option {letter}.")
            print("<<<{}>>>".format(letter))
            break