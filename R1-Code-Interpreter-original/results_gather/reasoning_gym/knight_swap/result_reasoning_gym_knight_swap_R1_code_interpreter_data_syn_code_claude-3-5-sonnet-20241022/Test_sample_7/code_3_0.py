from collections import deque
import copy

class State:
    def __init__(self):
        self.board = [['.'] * 3 for _ in range(4)]  # 4x3 board
        self.white_pos = set()
        self.black_pos = set()
        
    def initialize(self):
        # Set initial positions
        self.board[1][2] = 'B'  # B3
        self.board[2][1] = 'w'  # C2
        self.board[3][1] = 'w'  # D2
        self.board[2][0] = 'B'  # C1
        self.white_pos = {(2,1), (3,1)}
        self.black_pos = {(1,2), (2,0)}
        
    def is_valid(self, x, y):
        return 0 <= x < 4 and 0 <= y < 3
        
    def get_knight_moves(self, x, y):
        moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
                 (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
        return [(nx,ny) for nx,ny in moves if self.is_valid(nx,ny) and self.board[nx][ny] == '.']
        
    def make_move(self, from_pos, to_pos, is_white):
        x1,y1 = from_pos
        x2,y2 = to_pos
        piece = 'w' if is_white else 'B'
        self.board[x1][y1] = '.'
        self.board[x2][y2] = piece
        pos_set = self.white_pos if is_white else self.black_pos
        pos_set.remove((x1,y1))
        pos_set.add((x2,y2))
        
    def is_target_reached(self):
        return self.white_pos == {(1,2), (2,0)}  # Original black positions
        
    def get_state_key(self):
        return (tuple(sorted(self.white_pos)), tuple(sorted(self.black_pos)))

def coord_to_notation(x, y):
    return f"{chr(65+x)}{y+1}"

def solve():
    initial_state = State()
    initial_state.initialize()
    
    visited = set()
    queue = deque([(initial_state, [], True, 0)])  # state, moves, is_white_turn, depth
    best_solution = None
    max_depth = 7  # Limit search depth
    
    while queue:
        state, moves, is_white_turn, depth = queue.popleft()
        
        if depth >= max_depth:
            continue
            
        state_key = state.get_state_key()
        if state_key in visited:
            continue
        visited.add(state_key)
        
        if state.is_target_reached():
            if best_solution is None or len(moves) < len(best_solution):
                best_solution = moves
            continue
            
        current_positions = state.white_pos if is_white_turn else state.black_pos
        for x, y in current_positions:
            for new_x, new_y in state.get_knight_moves(x, y):
                new_state = copy.deepcopy(state)
                new_state.make_move((x,y), (new_x,new_y), is_white_turn)
                move = f"{'w' if is_white_turn else 'b'},{coord_to_notation(x,y)},{coord_to_notation(new_x,new_y)}"
                queue.append((new_state, moves + [move], not is_white_turn, depth + 1))
    
    return best_solution

# Solve and print result
result = solve()
if result:
    print(result)
else:
    print("No")