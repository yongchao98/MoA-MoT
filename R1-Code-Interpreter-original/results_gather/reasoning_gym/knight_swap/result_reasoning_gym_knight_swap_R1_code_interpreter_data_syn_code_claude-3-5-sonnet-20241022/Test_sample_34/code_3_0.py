from collections import deque
import copy

class BoardState:
    def __init__(self):
        # Initialize board with exact pieces
        self.board = [
            ['B', '.', '.', '.'],  # row 3 (index 0)
            ['B', '.', '.', 'w'],  # row 2 (index 1)
            ['.', 'w', '.', '.']   # row 1 (index 2)
        ]
        # Keep track of piece positions
        self.white_positions = {(1,2), (3,1)}  # B1, D2
        self.black_positions = {(0,0), (0,1)}  # A3, A2
        
    def clone(self):
        new_state = BoardState()
        new_state.board = copy.deepcopy(self.board)
        new_state.white_positions = self.white_positions.copy()
        new_state.black_positions = self.black_positions.copy()
        return new_state
    
    def make_move(self, from_pos, to_pos, is_white):
        x1, y1 = from_pos
        x2, y2 = to_pos
        
        # Update board
        self.board[y2][x2] = self.board[y1][x1]
        self.board[y1][x1] = '.'
        
        # Update position sets
        pos_set = self.white_positions if is_white else self.black_positions
        pos_set.remove((x1,y1))
        pos_set.add((x2,y2))
    
    def is_valid_move(self, from_pos, to_pos, is_white):
        x1, y1 = from_pos
        x2, y2 = to_pos
        
        # Check boundaries
        if not (0 <= x1 < 4 and 0 <= y1 < 3 and 0 <= x2 < 4 and 0 <= y2 < 3):
            return False
            
        # Check if moving correct piece
        if is_white and (x1,y1) not in self.white_positions:
            return False
        if not is_white and (x1,y1) not in self.black_positions:
            return False
            
        # Check if target square is empty
        if self.board[y2][x2] != '.':
            return False
            
        # Check L-shape move
        dx = abs(x2 - x1)
        dy = abs(y2 - y1)
        return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)
    
    def is_solved(self):
        # White should be at original black positions
        target_white = {(0,0), (0,1)}  # A3, A2
        # Black should be at original white positions
        target_black = {(1,2), (3,1)}  # B1, D2
        return self.white_positions == target_white and self.black_positions == target_black
    
    def get_state_key(self):
        return str(sorted(self.white_positions)) + str(sorted(self.black_positions))

def pos_to_notation(pos):
    return chr(pos[0] + ord('A')) + str(3 - pos[1])

def solve_knights():
    initial = BoardState()
    visited = set()
    queue = deque([(initial, [], True)])  # (state, moves, is_white_turn)
    
    while queue:
        current, moves, is_white = queue.popleft()
        state_key = current.get_state_key()
        
        if current.is_solved():
            return moves
            
        if state_key in visited:
            continue
            
        visited.add(state_key)
        
        # Get current pieces positions
        pieces = current.white_positions if is_white else current.black_positions
        
        # Try all possible moves for each piece
        for from_pos in pieces:
            x, y = from_pos
            possible_moves = [
                (x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
                (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)
            ]
            
            for to_pos in possible_moves:
                if current.is_valid_move(from_pos, to_pos, is_white):
                    # Create new state and make move
                    new_state = current.clone()
                    new_state.make_move(from_pos, to_pos, is_white)
                    
                    # Record move
                    piece = 'w' if is_white else 'B'
                    move = f"{piece},{pos_to_notation(from_pos)},{pos_to_notation(to_pos)}"
                    queue.append((new_state, moves + [move], not is_white))
    
    return "No"

result = solve_knights()
print(result)