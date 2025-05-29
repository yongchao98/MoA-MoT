from collections import deque
import copy

# Board representation
class Board:
    def __init__(self):
        self.size = 4
        self.board = [['' for _ in range(4)] for _ in range(3)]
        # Initial setup
        self.board[0][1] = 'w'
        self.board[0][2] = 'w'
        self.board[0][3] = 'B'
        self.board[1][0] = '.'
        self.board[1][3] = '.'
        self.board[2][1] = '.'
        self.board[2][2] = 'B'

    def __str__(self):
        return str(self.board)

    def get_state(self):
        return tuple(tuple(row) for row in self.board)

def is_valid_move(x1, y1, x2, y2, size):
    if not (0 <= x2 < 3 and 0 <= y2 < size):
        return False
    dx = abs(x2 - x1)
    dy = abs(y2 - y1)
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def get_valid_moves(board, is_white_turn):
    moves = []
    piece = 'w' if is_white_turn else 'B'
    for i in range(3):
        for j in range(4):
            if board.board[i][j] == piece:
                for x in range(3):
                    for y in range(4):
                        if (board.board[x][y] == '.' and 
                            is_valid_move(i, j, x, y, 4)):
                            moves.append((i, j, x, y))
    return moves

def is_goal_state(board):
    # Check if white knights are where black knights were and vice versa
    if (board.board[0][1] == 'B' and 
        board.board[0][2] == 'B' and 
        board.board[0][3] == 'w' and 
        board.board[2][2] == 'w'):
        return True
    return False

def solve():
    initial_board = Board()
    queue = deque([(initial_board, [], True)])  # (board, moves, is_white_turn)
    visited = set()
    visited.add(initial_board.get_state())

    while queue:
        current_board, moves, is_white_turn = queue.popleft()
        
        if is_goal_state(current_board):
            return moves
            
        valid_moves = get_valid_moves(current_board, is_white_turn)
        
        for move in valid_moves:
            new_board = copy.deepcopy(current_board)
            x1, y1, x2, y2 = move
            piece = new_board.board[x1][y1]
            new_board.board[x1][y1] = '.'
            new_board.board[x2][y2] = piece
            
            state = new_board.get_state()
            if state not in visited:
                visited.add(state)
                move_notation = f"{piece},{chr(65+y1)}{3-x1},{chr(65+y2)}{3-x2}"
                queue.append((new_board, moves + [move_notation], not is_white_turn))

    return None

# Solve the puzzle
solution = solve()
if solution:
    print(solution)
else:
    print("No")