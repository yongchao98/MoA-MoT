import copy

class RubiksCube:
    """
    A class to represent and manipulate a 3x3 Rubik's Cube.
    The cube state is stored in a dictionary of faces.
    F, B, U, D, L, R correspond to Front, Back, Up, Down, Left, Right faces.
    """
    def __init__(self, faces):
        self.faces = copy.deepcopy(faces)

    def _rotate_face_cw(self, face):
        return [[face[2][0], face[1][0], face[0][0]],
                [face[2][1], face[1][1], face[0][1]],
                [face[2][2], face[1][2], face[0][2]]]

    def _rotate_face_ccw(self, face):
        return [[face[0][2], face[1][2], face[2][2]],
                [face[0][1], face[1][1], face[2][1]],
                [face[0][0], face[1][0], face[2][0]]]

    def _move(self, face_key, is_clockwise):
        if is_clockwise:
            self.faces[face_key] = self._rotate_face_cw(self.faces[face_key])
        else:
            self.faces[face_key] = self._rotate_face_ccw(self.faces[face_key])

        # Store original strips before modification
        F_top, F_bottom = self.faces['F'][0][:], self.faces['F'][2][:]
        F_left = [row[0] for row in self.faces['F']]
        F_right = [row[2] for row in self.faces['F']]
        
        B_top, B_bottom = self.faces['B'][0][:], self.faces['B'][2][:]
        B_left = [row[0] for row in self.faces['B']]
        B_right = [row[2] for row in self.faces['B']]

        U_top, U_bottom = self.faces['U'][0][:], self.faces['U'][2][:]
        U_left = [row[0] for row in self.faces['U']]
        U_right = [row[2] for row in self.faces['U']]

        D_top, D_bottom = self.faces['D'][0][:], self.faces['D'][2][:]
        D_left = [row[0] for row in self.faces['D']]
        D_right = [row[2] for row in self.faces['D']]

        L_top, L_bottom = self.faces['L'][0][:], self.faces['L'][2][:]
        L_left = [row[0] for row in self.faces['L']]
        L_right = [row[2] for row in self.faces['L']]

        R_top, R_bottom = self.faces['R'][0][:], self.faces['R'][2][:]
        R_left = [row[0] for row in self.faces['R']]
        R_right = [row[2] for row in self.faces['R']]

        # CW: F->R->B->L->F
        # CCW: F->L->B->R->F
        if face_key == 'U':
            if is_clockwise:
                self.faces['F'][0], self.faces['R'][0], self.faces['B'][0], self.faces['L'][0] = L_top, F_top, R_top, B_top
            else: # U'
                self.faces['F'][0], self.faces['R'][0], self.faces['B'][0], self.faces['L'][0] = R_top, B_top, L_top, F_top
        
        # CW: F->L->B->R->F
        # CCW: F->R->B->L->F
        elif face_key == 'D':
            if is_clockwise:
                 self.faces['F'][2], self.faces['L'][2], self.faces['B'][2], self.faces['R'][2] = R_bottom, F_bottom, L_bottom, B_bottom
            else: # D'
                 self.faces['F'][2], self.faces['L'][2], self.faces['B'][2], self.faces['R'][2] = L_bottom, B_bottom, R_bottom, F_bottom

        # CW: U->F->D->B_rev->U
        # CCW: U->B_rev->D->F->U
        elif face_key == 'L':
            if is_clockwise:
                for i in range(3): self.faces['U'][i][0], self.faces['F'][i][0], self.faces['D'][i][0], self.faces['B'][2-i][2] = B_left[i], U_left[i], F_left[i], D_left[i]
            else: # L'
                for i in range(3): self.faces['U'][i][0], self.faces['F'][i][0], self.faces['D'][i][0], self.faces['B'][2-i][2] = F_left[i], D_left[i], B_left[i], U_left[i]

        # CW: U->B_rev->D->F->U
        # CCW: U->F->D->B_rev->U
        elif face_key == 'R':
            if is_clockwise:
                for i in range(3): self.faces['U'][i][2], self.faces['B'][2-i][0], self.faces['D'][i][2], self.faces['F'][i][2] = F_right[i], U_right[i], B_right[i], D_right[i]
            else: # R'
                for i in range(3): self.faces['U'][i][2], self.faces['B'][2-i][0], self.faces['D'][i][2], self.faces['F'][i][2] = B_right[i], D_right[i], F_right[i], U_right[i]
        
        # CW: U->R->D->L->U (with twists)
        # CCW: U->L->D->R->U (with twists)
        elif face_key == 'F':
            if is_clockwise:
                self.faces['U'][2] = [L_right[2], L_right[1], L_right[0]]
                for i in range(3): self.faces['R'][i][0] = U_bottom[i]
                self.faces['D'][0] = [R_left[2], R_left[1], R_left[0]]
                for i in range(3): self.faces['L'][i][2] = D_top[i]
            else: # F'
                self.faces['U'][2] = [R_left[0], R_left[1], R_left[2]]
                for i in range(3): self.faces['L'][i][2] = U_bottom[2-i]
                self.faces['D'][0] = [L_right[0], L_right[1], L_right[2]]
                for i in range(3): self.faces['R'][i][0] = D_top[2-i]
    
    def F(self): self._move('F', True)
    def F_prime(self): self._move('F', False)
    def U(self): self._move('U', True)
    def U_prime(self): self._move('U', False)
    def D(self): self._move('D', True)
    def D_prime(self): self._move('D', False)
    def R(self): self._move('R', True)
    def R_prime(self): self._move('R', False)
    def L(self): self._move('L', True)
    def L_prime(self): self._move('L', False)
    def B(self): self._move('B', True)
    def B_prime(self): self._move('B', False)

# --- Main execution ---
# Define faces based on standard orientation: F(white), U(orange), R(blue), L(green), D(red), B(yellow)
initial_layout = {
    'F': [['R','Y','Y'],['G','W','W'],['G','O','O']], # White
    'U': [['R','Y','W'],['B','O','Y'],['Y','R','O']], # Orange
    'R': [['G','G','G'],['R','B','B'],['B','B','O']], # Blue
    'B': [['Y','W','R'],['B','Y','O'],['G','O','B']], # Yellow
    'L': [['R','W','W'],['R','G','O'],['W','W','B']], # Green
    'D': [['B','Y','Y'],['R','R','G'],['W','G','O']], # Red
}

cube = RubiksCube(initial_layout)

# Based on the prompt's ambiguity, we execute R', U', F', L, D'
cube.R_prime()
cube.U_prime()
cube.F_prime()
cube.L()
cube.D_prime()

final_white_face = cube.faces['F']

print("The final state of the white face is:")
# Print the final matrix as requested for the user
print(final_white_face)
<<<E>>>