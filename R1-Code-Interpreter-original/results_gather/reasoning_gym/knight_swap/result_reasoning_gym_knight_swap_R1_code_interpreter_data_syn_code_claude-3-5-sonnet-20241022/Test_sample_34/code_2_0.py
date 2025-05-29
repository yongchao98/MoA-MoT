from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0<=x<4 and 0<=y<3]

def pos_to_notation(pos):
    return chr(pos[0] + ord('A')) + str(pos[1] + 1)

def notation_to_pos(notation):
    return (ord(notation[0]) - ord('A'), int(notation[1]) - 1)

class BoardState:
    def __init__(self):
        self.board = [
            ['B', '.', '.', '.'],  # row 3
            ['B', '.', '.', 'w'],  # row 2
            ['.', 'w', '.', '.']   # row 1
        ]
        
    def make_move(self, from_pos, to_pos):
        piece = self.board[from_pos[1]][from_pos[0]]
        self.board[from_pos[1]][from_pos[0]] = '.'
        self.board[to_pos[1]][to_pos[0]] = piece
        
    def is_valid_move(self, from_pos, to_pos, is_white_turn):
        if not (0 <= from_pos[0] < 4 and 0 <= from_pos[1] < 3):
            return False
        if not (0 <= to_pos[0] < 4 and 0 <= to_pos[1] < 3):
            return False
            
        piece = self.board[from_pos[1]][from_pos[0]]
        target = self.board[to_pos[1]][to_pos[0]]
        
        # Check if moving correct color piece to empty square
        if (is_white_turn and piece != 'w') or (not is_white_turn and piece != 'B'):
            return False
        if target != '.':
            return False
            
        # Check if move is L-shape
        dx = abs(to_pos[0] - from_pos[0])
        dy = abs(to_pos[1] - from_pos[1])
        return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)
        
    def is_target_reached(self):
        # White knights should be at A3 and A2
        # Black knights should be at B1 and D2
        return (self.board[2][0] == '.' and 
                self.board[2][1] == 'B' and 
                self.board[1][3] == 'B' and 
                self.board[2][2] == '.' and 
                self.board[2][3] == '.' and 
                self.board[0][0] == 'w' and 
                self.board[1][0] == 'w')
    
    def get_state_key(self):
        return ''.join(''.join(row) for row in self.board)
    
    def get_pieces(self, is_white):
        pieces = []
        piece = 'w' if is_white else 'B'
        for y in range(3):
            for x in range(4):
                if self.board[y][x] == piece:
                    pieces.append((x, y))
        return pieces

def solve_puzzle():
    initial_state = BoardState()
    visited = set()
    queue = deque([(initial_state, [], True)])  # state, moves, is_white_turn
    
    while queue:
        current_state, moves, is_white_turn = queue.popleft()
        state_key = current_state.get_state_key()
        
        if state_key in visited:
            continue
            
        visited.add(state_key)
        
        if current_state.is_target_reached():
            return moves
            
        # Get all pieces of current player
        pieces = current_state.get_pieces(is_white_turn)
        
        for from_pos in pieces:
            # Try all possible knight moves
            for to_pos in get_knight_moves(from_pos):
                if current_state.is_valid_move(from_pos, to_pos, is_white_turn):
                    # Create new state
                    new_state = BoardState()
                    new_state.board = copy.deepcopy(current_state.board)
                    new_state.make_move(from_pos, to_pos)
                    
                    # Record move
                    piece = 'w' if is_white_turn else 'B'
                    move = f"{piece},{pos_to_notation(from_pos)},{pos_to_notation(to_pos)}"
                    
                    queue.append((new_state, moves + [move], not is_white_turn))
    
    return "No"

result = solve_puzzle()
print(result)