import collections

class Cube:
    """
    Represents a Rubik's Cube and its transformations.
    Faces are indexed 0-5: Up, Right, Front, Down, Left, Back.
    Colors: W, R, G, Y, O, B (White, Red, Green, Yellow, Orange, Blue)
    Initial orientation: W on Up, G on Front, R on Right.
    """
    def __init__(self):
        self.faces = [
            ['W'] * 9,  # Up
            ['R'] * 9,  # Right
            ['G'] * 9,  # Front
            ['Y'] * 9,  # Down
            ['O'] * 9,  # Left
            ['B'] * 9   # Back
        ]

    def get_state(self):
        """Returns a hashable tuple representation of the cube state."""
        return tuple(tuple(face) for face in self.faces)

    def set_state(self, state):
        """Sets the cube state from a tuple."""
        self.faces = [list(face) for face in state]

    def _rotate_face_90(self, face_idx, clockwise=True):
        """Rotates a single face 90 degrees."""
        f = self.faces[face_idx]
        if clockwise:
            # (0,1,2,3,4,5,6,7,8) -> (6,3,0,7,4,1,8,5,2)
            self.faces[face_idx] = [f[6], f[3], f[0], f[7], f[4], f[1], f[8], f[5], f[2]]
        else:
            # (0,1,2,3,4,5,6,7,8) -> (2,5,8,1,4,7,0,3,6)
            self.faces[face_idx] = [f[2], f[5], f[8], f[1], f[4], f[7], f[0], f[3], f[6]]

    def move(self, move_char):
        """Applies a single move (e.g., 'R', 'U'', 'F2')."""
        moves = {
            'U': (0, [4, 2, 1, 5], [(0,1,2), (0,1,2), (0,1,2), (0,1,2)], True),
            'D': (3, [5, 1, 2, 4], [(6,7,8), (6,7,8), (6,7,8), (6,7,8)], True),
            'R': (1, [0, 2, 3, 5], [(2,5,8), (2,5,8), (2,5,8), (6,3,0)], True),
            'L': (4, [0, 5, 3, 2], [(0,3,6), (8,5,2), (0,3,6), (0,3,6)], True),
            'F': (2, [0, 4, 3, 1], [(6,7,8), (0,3,6), (2,1,0), (8,5,2)], True),
            'B': (5, [0, 1, 3, 4], [(0,1,2), (2,5,8), (8,7,6), (6,3,0)], True)
        }
        
        base_move = move_char[0]
        face_idx, adjacent_faces, slices, clockwise = moves[base_move]
        
        if len(move_char) == 1: # Clockwise
            rotations = 1
        elif move_char[1] == "'": # Counter-clockwise
            rotations = 3
            clockwise = not clockwise
        elif move_char[1] == '2': # 180 degrees
            rotations = 2
        
        for _ in range(rotations if clockwise else 1):
             self._rotate_face_90(face_idx, clockwise)
        
        for _ in range(rotations):
            f1_idx, f2_idx, f3_idx, f4_idx = adjacent_faces
            s1, s2, s3, s4 = slices
            
            f1, f2, f3, f4 = self.faces[f1_idx], self.faces[f2_idx], self.faces[f3_idx], self.faces[f4_idx]
            
            temp = [f1[i] for i in s1]
            for i_from, i_to in zip(s4, s1): f1[i_to] = f4[i_from]
            for i_from, i_to in zip(s3, s4): f4[i_to] = f3[i_from]
            for i_from, i_to in zip(s2, s3): f3[i_to] = f2[i_from]
            for i, val in enumerate(temp): f2[s2[i]] = val


    def apply_move_str(self, move_str):
        """Applies a sequence of moves from a string."""
        for move in move_str.split():
            self.move(move)

    def rotate_x(self): # Tip cube forward (like R L' M')
        # U->F, F->D, D->B(inv), B->U(inv)
        self._rotate_face_90(1, True) # R
        self._rotate_face_90(4, False) # L'
        old_U, old_F, old_D, old_B = self.faces[0], self.faces[2], self.faces[3], self.faces[5]
        self.faces[0] = old_F
        self.faces[2] = old_D
        self.faces[3] = old_B[::-1] # Inverted
        self.faces[5] = old_U[::-1] # Inverted

    def rotate_y(self): # Turn cube right (like U D' E')
        # F->L, L->B, B->R, R->F
        self._rotate_face_90(0, True) # U
        self._rotate_face_90(3, False) # D'
        old_R, old_F, old_L, old_B = self.faces[1], self.faces[2], self.faces[4], self.faces[5]
        self.faces[1] = old_B
        self.faces[2] = old_R
        self.faces[4] = old_F
        self.faces[5] = old_L


def count_solved_f2l_pairs(state):
    """
    Counts solved F2L pairs for a Y-top, O-front cube.
    U(Y), R(G), F(O), D(W), L(B), B(R)
    """
    faces = state
    U, R, F, D, L, B = faces[0], faces[1], faces[2], faces[3], faces[4], faces[5]
    
    # Centers
    W, G, O, Y, B, R_c = D[4], R[4], F[4], U[4], L[4], B[4]
    
    solved_count = 0
    
    # Slot FR: Corner(W,O,G), Edge(O,G)
    if D[2]==W and F[8]==O and R[6]==G and F[5]==O and R[3]==G:
        solved_count += 1
    # Slot FL: Corner(W,O,B), Edge(O,B)
    if D[0]==W and F[6]==O and L[8]==B and F[3]==O and L[5]==B:
        solved_count += 1
    # Slot BR: Corner(W,R,G), Edge(R,G)
    if D[8]==W and B[2]==R_c and R[8]==G and B[5]==R_c and R[5]==G:
        solved_count += 1
    # Slot BL: Corner(W,R,B), Edge(R,B)
    if D[6]==W and B[0]==R_c and L[2]==B and B[3]==R_c and L[3]==B:
        solved_count += 1
        
    return solved_count

def solve_f2l():
    """Main function to perform the scramble and search."""
    cube = Cube()
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D"
    cube.apply_move_str(scramble)

    # Reorient cube: x2 y (W-top, G-front -> Y-top, O-front)
    cube.rotate_x()
    cube.rotate_x()
    cube.rotate_y()
    
    initial_state = cube.get_state()
    
    if count_solved_f2l_pairs(initial_state) >= 2:
        print(0)
        return

    q = collections.deque([(initial_state, 0, None)]) # state, depth, last_move_axis
    visited = {initial_state}
    
    moves = ["U", "U'", "U2", "D", "D'", "D2", "R", "R'", "R2", 
             "L", "L'", "L2", "F", "F'", "F2", "B", "B'", "B2"]
    
    move_axis = {m[0] for m in moves}
    
    while q:
        current_state, depth, last_axis = q.popleft()
        
        # Heuristic depth limit to prevent infinite search on error
        if depth >= 15:
            continue
            
        for move in moves:
            # Optimization: don't turn the same face twice in a row (e.g. R R')
            current_axis = move[0]
            if current_axis == last_axis:
                continue

            temp_cube = Cube()
            temp_cube.set_state(current_state)
            temp_cube.move(move)
            next_state = temp_cube.get_state()

            if next_state not in visited:
                if count_solved_f2l_pairs(next_state) >= 2:
                    print(depth + 1)
                    return
                
                visited.add(next_state)
                q.append((next_state, depth + 1, current_axis))

if __name__ == '__main__':
    solve_f2l()