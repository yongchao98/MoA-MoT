from collections import deque
import copy

# Board representation
class Board:
    def __init__(self):
        self.size = 4
        self.board = [['.'] * 4 for _ in range(4)]
        # Initial position
        self.board[2][1] = 'w'  # A2
        self.board[2][0] = 'w'  # B2
        self.board[0][1] = 'w'  # B3
        self.board[0][3] = 'B'  # D3
        self.board[0][2] = '.'  # C3
        self.board[1][1] = '.'  # B2
        self.board[1][2] = '.'  # C2
        self.board[1][3] = '.'  # D2
        self.board[2][2] = 'B'  # C1

    def is_valid_pos(self, x, y):
        return 0 <= x < self.size and 0 <= y < self.size

    def get_knight_moves(self, x, y):
        moves = [
            (x+2, y+1), (x+2, y-1),
            (x-2, y+1), (x-2, y-1),
            (x+1, y+2), (x+1, y-2),
            (x-1, y+2), (x-1, y-2)
        ]
        return [(nx, ny) for nx, ny in moves if self.is_valid_pos(nx, ny)]

    def get_state_key(self):
        return ''.join(''.join(row) for row in self.board)

    def is_goal_state(self):
        # Check if white knights are where black knights were and vice versa
        white_positions = {(2,1), (2,0), (0,1)}
        black_positions = {(0,3), (2,2)}
        
        for i in range(self.size):
            for j in range(self.size):
                if self.board[i][j] == 'w' and (i,j) not in black_positions:
                    return False
                if self.board[i][j] == 'B' and (i,j) not in white_positions:
                    return False
        return True

def solve_knight_swap():
    initial_board = Board()
    queue = deque([(initial_board, [], 'w')])  # board, moves, next_turn
    visited = {initial_board.get_state_key()}

    while queue:
        current_board, moves, turn = queue.popleft()
        
        if current_board.is_goal_state():
            return moves

        # Find all pieces of current color
        for i in range(current_board.size):
            for j in range(current_board.size):
                if current_board.board[i][j] == turn:
                    # Get all possible moves for this piece
                    for new_x, new_y in current_board.get_knight_moves(i, j):
                        if current_board.board[new_x][new_y] == '.':
                            # Make move
                            new_board = copy.deepcopy(current_board)
                            new_board.board[new_x][new_y] = turn
                            new_board.board[i][j] = '.'
                            
                            new_state = new_board.get_state_key()
                            if new_state not in visited:
                                visited.add(new_state)
                                # Convert coordinates to chess notation
                                from_pos = f"{chr(j+65)}{3-i}"
                                to_pos = f"{chr(new_y+65)}{3-new_x}"
                                new_moves = moves + [f"{turn},{from_pos},{to_pos}"]
                                queue.append((new_board, new_moves, 'B' if turn == 'w' else 'w'))

    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")