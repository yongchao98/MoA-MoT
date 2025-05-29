from collections import deque

class BoardState:
    def __init__(self):
        # Board representation: 0=empty, 1=white knight, 2=black knight
        self.board = [
            [2, 0, 0, 0],  # row 3 (A3,B3,C3,D3)
            [2, 0, 0, 1],  # row 2 (A2,B2,C2,D2)
            [0, 1, 0, 0]   # row 1 (A1,B1,C1,D1)
        ]
        
    def clone(self):
        new_state = BoardState()
        new_state.board = [row[:] for row in self.board]
        return new_state
    
    def is_valid_move(self, x1, y1, x2, y2, is_white):
        # Check boundaries
        if not (0 <= x1 < 4 and 0 <= y1 < 3 and 0 <= x2 < 4 and 0 <= y2 < 3):
            return False
        
        # Check if moving correct piece and to empty square
        piece = 1 if is_white else 2
        if self.board[y1][x1] != piece or self.board[y2][x2] != 0:
            return False
        
        # Check L-shape move
        dx = abs(x2 - x1)
        dy = abs(y2 - y1)
        return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)
    
    def make_move(self, x1, y1, x2, y2):
        piece = self.board[y1][x1]
        self.board[y1][x1] = 0
        self.board[y2][x2] = piece
    
    def is_solved(self):
        # White knights should be at A3,A2; Black knights at B1,D2
        return (self.board[0][0] == 1 and  # A3
                self.board[1][0] == 1 and  # A2
                self.board[2][1] == 2 and  # B1
                self.board[1][3] == 2)     # D2
    
    def get_state_key(self):
        return str(self.board)

def get_knight_moves(x, y):
    moves = [
        (x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
        (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)
    ]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def pos_to_notation(x, y):
    return f"{chr(x + ord('A'))}{3-y}"

def solve():
    start = BoardState()
    visited = set()
    queue = deque([(start, [], True)])  # (state, moves, is_white_turn)
    
    while queue:
        current, moves, is_white = queue.popleft()
        
        if current.is_solved():
            return moves
        
        state_key = current.get_state_key()
        if state_key in visited:
            continue
        visited.add(state_key)
        
        piece = 1 if is_white else 2
        # Find all pieces of current color
        for y in range(3):
            for x in range(4):
                if current.board[y][x] == piece:
                    # Try all possible moves
                    for new_x, new_y in get_knight_moves(x, y):
                        if current.is_valid_move(x, y, new_x, new_y, is_white):
                            new_state = current.clone()
                            new_state.make_move(x, y, new_x, new_y)
                            
                            move = f"{'w' if is_white else 'B'},{pos_to_notation(x,y)},{pos_to_notation(new_x,new_y)}"
                            if len(moves) < 6:  # Limit depth to prevent long solutions
                                queue.append((new_state, moves + [move], not is_white))
    
    return "No"

result = solve()
print(result)