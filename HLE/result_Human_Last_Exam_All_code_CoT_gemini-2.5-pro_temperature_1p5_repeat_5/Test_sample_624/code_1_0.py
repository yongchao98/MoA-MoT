import collections
import copy

class RubiksCube:
    """
    A class to represent and manipulate a Rubik's Cube.
    The state is stored as a dictionary of faces, with standard colors:
    U: White (W), F: Green (G), R: Red (R), L: Orange (O), B: Blue (B), D: Yellow (Y)
    """

    def __init__(self, state_tuple=None):
        if state_tuple is None:
            self.state = {
                'U': [['W'] * 3 for _ in range(3)], 'D': [['Y'] * 3 for _ in range(3)],
                'F': [['G'] * 3 for _ in range(3)], 'B': [['B'] * 3 for _ in range(3)],
                'L': [['O'] * 3 for _ in range(3)], 'R': [['R'] * 3 for _ in range(3)]
            }
        else:
            self.state = {
                'U': [list(r) for r in state_tuple[0]], 'L': [list(r) for r in state_tuple[1]],
                'F': [list(r) for r in state_tuple[2]], 'R': [list(r) for r in state_tuple[3]],
                'B': [list(r) for r in state_tuple[4]], 'D': [list(r) for r in state_tuple[5]]
            }

    @staticmethod
    def _rotate_face(face, clockwise=True):
        if clockwise:
            return [list(row) for row in zip(*face[::-1])]
        else: # counter-clockwise
            return [list(row) for row in zip(*face)][::-1]

    def to_tuple(self):
        return tuple(tuple(map(tuple, self.state[f])) for f in ['U', 'L', 'F', 'R', 'B', 'D'])

    def U(self, n=1):
        for _ in range(n % 4):
            self.state['U'] = self._rotate_face(self.state['U'])
            temp_row = self.state['F'][0]
            self.state['F'][0] = self.state['R'][0]
            self.state['R'][0] = self.state['B'][0]
            self.state['B'][0] = self.state['L'][0]
            self.state['L'][0] = temp_row

    def D(self, n=1):
        for _ in range(n % 4):
            self.state['D'] = self._rotate_face(self.state['D'])
            temp_row = self.state['F'][2]
            self.state['F'][2] = self.state['L'][2]
            self.state['L'][2] = self.state['B'][2]
            self.state['B'][2] = self.state['R'][2]
            self.state['R'][2] = temp_row

    def R(self, n=1):
        for _ in range(n % 4):
            self.state['R'] = self._rotate_face(self.state['R'])
            temp_col = [row[2] for row in self.state['U']]
            for i in range(3): self.state['U'][i][2] = self.state['F'][i][2]
            for i in range(3): self.state['F'][i][2] = self.state['D'][i][2]
            for i in range(3): self.state['D'][i][2] = self.state['B'][2 - i][0]
            for i in range(3): self.state['B'][2 - i][0] = temp_col[i]
            
    def L(self, n=1):
        for _ in range(n % 4):
            self.state['L'] = self._rotate_face(self.state['L'])
            temp_col = [row[0] for row in self.state['U']]
            for i in range(3): self.state['U'][i][0] = self.state['B'][2 - i][2]
            for i in range(3): self.state['B'][2-i][2] = self.state['D'][i][0]
            for i in range(3): self.state['D'][i][0] = self.state['F'][i][0]
            for i in range(3): self.state['F'][i][0] = temp_col[i]

    def F(self, n=1):
        for _ in range(n % 4):
            self.state['F'] = self._rotate_face(self.state['F'])
            temp_row = self.state['U'][2]
            self.state['U'][2] = [self.state['L'][2 - i][2] for i in range(3)]
            for i in range(3): self.state['L'][i][2] = self.state['D'][0][i]
            self.state['D'][0] = [self.state['R'][2 - i][0] for i in range(3)]
            for i in range(3): self.state['R'][i][0] = temp_row[i]

    def B(self, n=1):
        for _ in range(n % 4):
            self.state['B'] = self._rotate_face(self.state['B'])
            temp_row = self.state['U'][0]
            self.state['U'][0] = [self.state['R'][i][2] for i in range(3)]
            for i in range(3): self.state['R'][i][2] = self.state['D'][2][2-i]
            self.state['D'][2] = [self.state['L'][i][0] for i in range(3)]
            for i in range(3): self.state['L'][i][0] = temp_row[2-i]

    def apply_move(self, move):
        n = 1
        if len(move) > 1:
            if move[1] == "'": n = 3
            elif move[1] == "2": n = 2
        
        move_func = getattr(self, move[0])
        move_func(n)

def is_goal(cube):
    solved_pairs = 0
    state = cube.state
    
    # New orientation (Y-top, O-front) corresponds to centers:
    # U=Y, D=W, F=O, R=B, B=R, L=G
    # Standard cube colors are: W,G,R,B,O,Y
    
    # Check FR slot (Corner: WOB, Edge: OB)
    if (state['D'][0][2] == 'W' and state['F'][2][2] == 'O' and state['R'][2][0] == 'B' and
        state['F'][1][2] == 'O' and state['R'][1][0] == 'B'):
        solved_pairs += 1
        
    # Check FL slot (Corner: WOG, Edge: OG)
    if (state['D'][0][0] == 'W' and state['F'][2][0] == 'O' and state['L'][2][2] == 'G' and
        state['F'][1][0] == 'O' and state['L'][1][2] == 'G'):
        solved_pairs += 1
        
    # Check BR slot (Corner: WRB, Edge: RB)
    if (state['D'][2][2] == 'W' and state['B'][2][0] == 'R' and state['R'][2][2] == 'B' and
        state['B'][1][0] == 'R' and state['R'][1][2] == 'B'):
        solved_pairs += 1

    # Check BL slot (Corner: WRG, Edge: RG)
    if (state['D'][2][0] == 'W' and state['B'][2][2] == 'R' and state['L'][2][0] == 'G' and
        state['B'][1][2] == 'R' and state['L'][1][0] == 'G'):
        solved_pairs += 1
        
    return solved_pairs >= 2

def solve_two_f2l(start_cube):
    # Moves from new orientation (Y-top, O-front) mapped to original (W-top, G-front)
    # New U -> Old D, New F -> Old L, New R -> Old F, etc. (rotation x2 y)
    move_map = {'U': 'D', 'D': 'U', 'F': 'L', 'L': 'B', 'B': 'R', 'R': 'F'}
    
    # All possible moves in the new orientation
    moves_to_try = []
    for base in move_map.keys():
        moves_to_try.extend([base, base + "'", base + "2"])
        
    queue = collections.deque([(start_cube.to_tuple(), [])])
    visited = {start_cube.to_tuple()}
    
    while queue:
        current_state_tuple, path = queue.popleft()
        current_cube = RubiksCube(current_state_tuple)

        if is_goal(current_cube):
            return path
            
        for move_str in moves_to_try:
            next_cube = RubiksCube(current_state_tuple)
            
            # Apply the move based on the mapping
            original_move_char = move_map[move_str[0]]
            move_to_apply = original_move_char + (move_str[1:] if len(move_str) > 1 else "")
            next_cube.apply_move(move_to_apply)

            next_state_tuple = next_cube.to_tuple()
            if next_state_tuple not in visited:
                visited.add(next_state_tuple)
                new_path = path + [move_str]
                queue.append((next_state_tuple, new_path))
    return None

if __name__ == '__main__':
    scramble = "R2 U F2 U2 B2 U2 R' U2 F2 R' B2 D2 R U2 D L' R2 D2 F L B' F' U' R L' D R' B R F2 L' F2 L D"
    
    cube = RubiksCube()
    for move in scramble.split():
        cube.apply_move(move)

    solution_path = solve_two_f2l(cube)

    if solution_path:
        # Each move is printed as part of the equation
        print(" ".join(solution_path))
        print(f"<<<{len(solution_path)}>>>")
    else:
        print("No solution found.")