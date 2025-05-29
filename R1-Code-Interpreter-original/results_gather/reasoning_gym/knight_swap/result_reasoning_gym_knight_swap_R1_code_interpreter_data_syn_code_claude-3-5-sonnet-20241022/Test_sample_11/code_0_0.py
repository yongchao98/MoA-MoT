from collections import deque
import copy

# Board representation
class Board:
    def __init__(self):
        self.size = 4
        self.board = [['' for _ in range(4)] for _ in range(4)]
        # Initial setup
        self.board[2][1] = 'w'  # B3
        self.board[1][0] = 'w'  # A2
        self.board[1][3] = 'B'  # D2
        self.board[0][3] = 'B'  # D1

    def is_valid_pos(self, x, y):
        return 0 <= x < self.size and 0 <= y < self.size

    def get_knight_moves(self, x, y):
        moves = [
            (x+2, y+1), (x+2, y-1),
            (x-2, y+1), (x-2, y-1),
            (x+1, y+2), (x+1, y-2),
            (x-1, y+2), (x-1, y-2)
        ]
        return [(nx, ny) for nx, ny in moves if self.is_valid_pos(nx, ny) and self.board[nx][ny] == '']

    def to_string(self):
        return ''.join(''.join(row) for row in self.board)

    def is_goal_state(self):
        # Check if white knights are where black knights started and vice versa
        black_positions = {(1,3), (0,3)}
        white_positions = {(2,1), (1,0)}
        
        for i in range(self.size):
            for j in range(self.size):
                if self.board[i][j] == 'w' and (i,j) not in black_positions:
                    return False
                if self.board[i][j] == 'B' and (i,j) not in white_positions:
                    return False
        return True

def solve_knight_swap():
    initial_board = Board()
    queue = deque([(initial_board, [], 'w')])  # (board, moves, next_turn)
    visited = {initial_board.to_string()}

    while queue:
        current_board, moves, turn = queue.popleft()
        
        if current_board.is_goal_state():
            return moves

        # Find all knights of current color
        knights = []
        for i in range(current_board.size):
            for j in range(current_board.size):
                if current_board.board[i][j] == turn:
                    knights.append((i, j))

        # Try all possible moves for each knight
        for kx, ky in knights:
            valid_moves = current_board.get_knight_moves(kx, ky)
            for nx, ny in valid_moves:
                new_board = Board()
                new_board.board = copy.deepcopy(current_board.board)
                new_board.board[nx][ny] = new_board.board[kx][ky]
                new_board.board[kx][ky] = ''
                
                board_str = new_board.to_string()
                if board_str not in visited:
                    visited.add(board_str)
                    # Convert coordinates to chess notation
                    from_pos = f"{chr(ky+65)}{4-kx}"
                    to_pos = f"{chr(ny+65)}{4-nx}"
                    new_moves = moves + [f"{turn},{from_pos},{to_pos}"]
                    queue.append((new_board, new_moves, 'B' if turn == 'w' else 'w'))

    return None

# Solve and print result
result = solve_knight_swap()
if result is None:
    print("<<<No>>>")
else:
    print(f"<<<{result}>>>")