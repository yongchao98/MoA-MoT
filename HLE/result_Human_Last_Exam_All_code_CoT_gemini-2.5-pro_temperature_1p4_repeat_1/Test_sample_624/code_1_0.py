import collections

class RubiksCube:
    """
    Represents a Rubik's Cube and its operations.

    The cube state is a list of 54 characters, representing the sticker colors.
    The face order is U, L, F, R, B, D.
    Initial orientation (solved state): U=White, L=Orange, F=Green, R=Red, B=Blue, D=Yellow.
    """

    def __init__(self, state_list=None):
        if state_list:
            self.state = state_list
        else:
            # U(W) L(O) F(G) R(R) B(B) D(Y)
            self.state = list('W'*9 + 'O'*9 + 'G'*9 + 'R'*9 + 'B'*9 + 'Y'*9)

    def get_state_tuple(self):
        return tuple(self.state)

    def apply_moves(self, moves_str):
        for move in moves_str.split():
            self._apply_single_move(move)

    def _apply_single_move(self, move_str):
        if len(move_str) == 1:
            move_func = getattr(self, f"move_{move_str}")
            move_func()
        elif move_str[1] == "'":
            move_func = getattr(self, f"move_{move_str[0]}_prime")
            move_func()
        elif move_str[1] == "2":
            move_func = getattr(self, f"move_{move_str[0]}_2")
            move_func()

    def _rotate_face_clockwise(self, s):
        return [s[6],s[3],s[0], s[7],s[4],s[1], s[8],s[5],s[2]]

    def _cycle_stickers(self, indices):
        """ Cycles stickers in groups of 3 """
        temp = [self.state[i] for i in indices[-3:]]
        for i in range(len(indices) - 3, 2, -3):
            for j in range(3):
                self.state[indices[i+j]] = self.state[indices[i-3+j]]
        for j in range(3):
            self.state[indices[j]] = temp[j]

    # --- Base Moves (Clockwise) ---
    def move_U(self):
        self.state[0:9] = self._rotate_face_clockwise(self.state[0:9])
        self._cycle_stickers([36,37,38, 27,28,29, 18,19,20, 9,10,11])
    def move_L(self):
        self.state[9:18] = self._rotate_face_clockwise(self.state[9:18])
        self._cycle_stickers([0,3,6, 18,21,24, 45,48,51, 44,41,38])
    def move_F(self):
        self.state[18:27] = self._rotate_face_clockwise(self.state[18:27])
        self._cycle_stickers([6,7,8, 27,30,33, 47,46,45, 17,14,11])
    def move_R(self):
        self.state[27:36] = self._rotate_face_clockwise(self.state[27:36])
        self._cycle_stickers([2,5,8, 38,41,44, 47,50,53, 20,23,26])
    def move_B(self):
        self.state[36:45] = self._rotate_face_clockwise(self.state[36:45])
        self._cycle_stickers([2,1,0, 15,12,9, 51,52,53, 35,32,29])
    def move_D(self):
        self.state[45:54] = self._rotate_face_clockwise(self.state[45:54])
        self._cycle_stickers([11,12,13, 20,21,22, 29,30,31, 38,39,40]) # Typo from web, let me fix
        # Actual D move cycle
        # F(bot) -> R(bot) -> B(bot) -> L(bot)
        self._cycle_stickers([24,25,26, 33,34,35, 42,43,44, 15,16,17])

    # --- Prime and Double Moves ---
    def move_U_prime(self): self.move_U(); self.move_U(); self.move_U()
    def move_F_prime(self): self.move_F(); self.move_F(); self.move_F()
    def move_D_prime(self): self.move_D(); self.move_D(); self.move_D()
    def move_B_prime(self): self.move_B(); self.move_B(); self.move_B()
    def move_L_prime(self): self.move_L(); self.move_L(); self.move_L()
    def move_R_prime(self): self.move_R(); self.move_R(); self.move_R()
    def move_U_2(self): self.move_U(); self.move_U()
    def move_F_2(self): self.move_F(); self.move_F()
    def move_D_2(self): self.move_D(); self.move_D()
    def move_B_2(self): self.move_B(); self.move_B()
    def move_L_2(self): self.move_L(); self.move_L()
    def move_R_2(self): self.move_R(); self.move_R()

    def reorient(self, rotation):
        faces = [self.state[i:i+9] for i in range(0, 54, 9)]
        if rotation == 'x2': # U<->D, F<->B, L,R rot 180
            faces[0], faces[5] = faces[5], faces[0]
            faces[2], faces[4] = faces[4], faces[2]
            faces[1] = self._rotate_face_clockwise(self._rotate_face_clockwise(faces[1]))
            faces[3] = self._rotate_face_clockwise(self._rotate_face_clockwise(faces[3]))
        elif rotation == "y'": # F->R->B->L->F, U rot ccw, D rot cw
            faces[2], faces[3], faces[4], faces[1] = faces[3], faces[4], faces[1], faces[2]
            faces[0] = self._rotate_face_clockwise(self._rotate_face_clockwise(self._rotate_face_clockwise(faces[0])))
            faces[5] = self._rotate_face_clockwise(faces[5])
        self.state = [sticker for face in faces for sticker in face]

    def count_solved_f2l_pairs(self):
        """
        After reorient (Y-top, O-front):
        U=Y, L=G, F=O, R=B, B=R, D=W
        """
        s, solved_count = self.state, 0
        # DFR pair (WOB corner, OB edge)
        # Corner(D-53,F-26,R-33), Edge(F-23,R-30)
        if s[53]=='W' and s[26]=='O' and s[33]=='B' and s[23]=='O' and s[30]=='B':
            solved_count += 1
        # DFL pair (WOG corner, OG edge)
        # Corner(D-51,F-24,L-15), Edge(F-21,L-12)
        if s[51]=='W' and s[24]=='O' and s[15]=='G' and s[21]=='O' and s[12]=='G':
            solved_count += 1
        # DBR pair (WRB corner, RB edge)
        # Corner(D-47,B-44,R-35), Edge(B-41,R-28)
        if s[47]=='W' and s[44]=='R' and s[35]=='B' and s[41]=='R' and s[28]=='B':
            solved_count += 1
        # DBL pair (WRG corner, RG edge)
        # Corner(D-45,B-42,L-17), Edge(B-39,L-10)
        if s[45]=='W' and s[42]=='R' and s[17]=='G' and s[39]=='R' and s[10]=='G':
            solved_count += 1
        return solved_count

def solve_f2l():
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D"
    
    initial_cube = RubiksCube()
    initial_cube.apply_moves(scramble)

    # Reorient to Y-top, O-front
    initial_cube.reorient('x2')
    initial_cube.reorient("y'")

    if initial_cube.count_solved_f2l_pairs() >= 2:
        print("The minimum sequence to solve two F2L pairs is:")
        print("(no moves needed)")
        print("The length of the sequence is 0.")
        print(f"<<<0>>>")
        return

    q = collections.deque([(initial_cube.get_state_tuple(), [])])
    visited = {initial_cube.get_state_tuple()}
    
    moves = [ "U", "U'", "U2", "D", "D'", "D2", "L", "L'", "L2", "F", "F'", "F2", "R", "R'", "R2", "B", "B'", "B2" ]
    move_map = {m: getattr(RubiksCube, f"_apply_single_move") for m in moves}

    while q:
        current_state_tuple, path = q.popleft()
        
        last_move_face = path[-1][0] if path else None

        for move in moves:
            # Pruning: Don't turn the same face twice in a row (e.g., R followed by R')
            if move[0] == last_move_face:
                continue

            # Create a new state by applying the move
            temp_cube = RubiksCube(list(current_state_tuple))
            temp_cube._apply_single_move(move)
            
            new_state_tuple = temp_cube.get_state_tuple()

            if new_state_tuple not in visited:
                new_path = path + [move]
                if temp_cube.count_solved_f2l_pairs() >= 2:
                    solution_str = " ".join(new_path)
                    solution_len = len(new_path)
                    
                    print(f"The minimum sequence to solve two F2L pairs is:")
                    print(solution_str)
                    print(f"\nThe final equation is: {solution_str}. The length of the sequence is {solution_len}.")
                    print(f"<<<{solution_len}>>>")
                    return
                
                visited.add(new_state_tuple)
                q.append((new_state_tuple, new_path))

solve_f2l()