import copy

class RubiksCube:
    """
    A class to represent and manipulate a 3x3 Rubik's Cube.
    """

    def __init__(self):
        """
        Initializes the cube to the specified jumbled state.
        F(White), U(Orange), R(Blue), B(Yellow), L(Green), D(Red)
        """
        self.faces = {
            'F': [['R', 'Y', 'Y'], ['G', 'W', 'W'], ['G', 'O', 'O']], # Front (White)
            'U': [['R', 'Y', 'W'], ['B', 'O', 'Y'], ['Y', 'R', 'O']], # Up (Orange)
            'R': [['G', 'G', 'G'], ['R', 'B', 'B'], ['B', 'B', 'O']], # Right (Blue)
            'B': [['Y', 'W', 'R'], ['B', 'Y', 'O'], ['G', 'O', 'B']], # Back (Yellow)
            'L': [['R', 'W', 'W'], ['R', 'G', 'O'], ['W', 'W', 'B']], # Left (Green)
            'D': [['B', 'Y', 'Y'], ['R', 'R', 'G'], ['W', 'G', 'O']], # Down (Red)
        }

    def _rotate_face_clockwise(self, face):
        """Rotates a 3x3 face 90 degrees clockwise."""
        return [list(row) for row in zip(*face[::-1])]

    def _rotate_face_counter_clockwise(self, face):
        """Rotates a 3x3 face 90 degrees counter-clockwise."""
        return [list(row) for row in reversed(list(zip(*face)))]

    def R(self):
        """Performs a clockwise R (Right) move."""
        self.faces['R'] = self._rotate_face_clockwise(self.faces['R'])
        
        # In an R move, the cycle is F -> U -> B -> D -> F for the affected columns.
        temp_col = [row[2] for row in self.faces['F']] # Store F's right column
        
        # D -> F
        for i in range(3): self.faces['F'][i][2] = self.faces['D'][i][2]
        # B -> D (twisted)
        self.faces['D'][0][2] = self.faces['B'][2][0]
        self.faces['D'][1][2] = self.faces['B'][1][0]
        self.faces['D'][2][2] = self.faces['B'][0][0]
        # U -> B (twisted)
        self.faces['U'][0][2] = self.faces['B'][2][0] # This uses B before it's updated, so we need to store it
        temp_b_col = [self.faces['B'][i][0] for i in range(3)]
        self.faces['B'][0][0] = self.faces['U'][2][2]
        self.faces['B'][1][0] = self.faces['U'][1][2]
        self.faces['B'][2][0] = self.faces['U'][0][2]
        # F (temp) -> U
        for i in range(3): self.faces['U'][i][2] = temp_col[i]
        
    def U(self):
        """Performs a clockwise U (Up) move."""
        self.faces['U'] = self._rotate_face_clockwise(self.faces['U'])
        # Cycle is F -> R -> B -> L -> F for top rows
        temp_row = copy.deepcopy(self.faces['F'][0])
        self.faces['F'][0] = self.faces['R'][0]
        self.faces['R'][0] = self.faces['B'][0]
        self.faces['B'][0] = self.faces['L'][0]
        self.faces['L'][0] = temp_row

    def F(self):
        """Performs a clockwise F (Front) move."""
        self.faces['F'] = self._rotate_face_clockwise(self.faces['F'])
        # Cycle is U -> R -> D -> L -> U
        temp_row = copy.deepcopy(self.faces['U'][2])
        # L right col -> U bottom row (twisted)
        self.faces['U'][2] = [self.faces['L'][2][2], self.faces['L'][1][2], self.faces['L'][0][2]]
        # D top row -> L right col
        self.faces['L'][0][2], self.faces['L'][1][2], self.faces['L'][2][2] = self.faces['D'][0]
        # R left col -> D top row (twisted)
        self.faces['D'][0] = [self.faces['R'][2][0], self.faces['R'][1][0], self.faces['R'][0][0]]
        # U bottom row (temp) -> R left col
        self.faces['R'][0][0], self.faces['R'][1][0], self.faces['R'][2][0] = temp_row

    def L_prime(self):
        """Performs a counter-clockwise L' (Left) move."""
        self.faces['L'] = self._rotate_face_counter_clockwise(self.faces['L'])
        # Cycle is U <- F <- D <- B <- U
        temp_col = [row[0] for row in self.faces['U']] # store U's left column
        # F -> U
        for i in range(3): self.faces['U'][i][0] = self.faces['F'][i][0]
        # D -> F
        for i in range(3): self.faces['F'][i][0] = self.faces['D'][i][0]
        # B -> D (twisted)
        self.faces['D'][0][0] = self.faces['B'][2][2]
        self.faces['D'][1][0] = self.faces['B'][1][2]
        self.faces['D'][2][0] = self.faces['B'][0][2]
        # U (temp) -> B (twisted)
        self.faces['B'][0][2] = temp_col[2]
        self.faces['B'][1][2] = temp_col[1]
        self.faces['B'][2][2] = temp_col[0]

    def D(self):
        """Performs a clockwise D (Down) move."""
        self.faces['D'] = self._rotate_face_clockwise(self.faces['D'])
        # Cycle is F -> L -> B -> R -> F for bottom rows
        temp_row = copy.deepcopy(self.faces['F'][2])
        self.faces['F'][2] = self.faces['L'][2]
        self.faces['L'][2] = self.faces['B'][2]
        self.faces['B'][2] = self.faces['R'][2]
        self.faces['R'][2] = temp_row

    def apply_moves(self, moves_str):
        """Applies a sequence of moves from a string."""
        moves = moves_str.split()
        for move in moves:
            if move == 'R': self.R()
            elif move == 'U': self.U()
            elif move == 'F': self.F()
            elif move == "L'": self.L_prime()
            elif move == 'D': self.D()

    def print_face(self, face_name):
        """Prints a single face of the cube."""
        print(f"Final state of the {face_name} face:")
        face = self.faces[face_name]
        # We need to print each character in the final list of lists
        print("[", end="")
        for i, row in enumerate(face):
            print(f"[{','.join(row)}]", end="")
            if i < len(face) - 1:
                print(",", end="")
        print("]")

if __name__ == '__main__':
    # 1. Create a cube instance
    cube = RubiksCube()

    # 2. Apply the algorithm R U F L' D
    cube.apply_moves("R U F L' D")

    # 3. Print the final state of the Front (White) face
    cube.print_face('F')