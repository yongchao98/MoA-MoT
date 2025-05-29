from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_valid_move(board, from_x, from_y, to_x, to_y, piece):
    # Check if source has the right piece
    if board[from_y][from_x] != piece:
        return False
    # Check if destination is empty
    if board[to_y][to_x] != '.':
        return False
    # Check if it's a valid knight move
    dx = abs(to_x - from_x)
    dy = abs(to_y - from_y)
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def is_goal_state(board):
    # White knights should be at B1, D1
    # Black knights should be at A2, C3
    return (board[2][1] == 'w' and board[2][3] == 'w' and 
            board[1][0] == 'B' and board[0][2] == 'B')

def find_solution():
    initial_board = [
        ['.', ' ', 'w', '.'],  # row 3
        ['w', ' ', '.', ' '],  # row 2
        ['.', 'B', '.', 'B']   # row 1
    ]
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    visited = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        piece = 'w' if is_white_turn else 'B'
        
        # Find all pieces of current color
        pieces = []
        for y in range(3):
            for x in range(4):
                if board[y][x] == piece:
                    pieces.append((x,y))
        
        # Try all possible moves for each piece
        for x,y in pieces:
            for new_x, new_y in get_knight_moves((x,y)):
                if is_valid_move(board, x, y, new_x, new_y, piece):
                    new_board = copy.deepcopy(board)
                    new_board[y][x] = '.'
                    new_board[new_y][new_x] = piece
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        from_pos = f"{chr(65+x)}{y+1}"
                        to_pos = f"{chr(65+new_x)}{new_y+1}"
                        new_moves = moves + [f"{piece},{from_pos},{to_pos}"]
                        
                        if is_goal_state(new_board):
                            return new_moves
                        
                        queue.append((new_board, new_moves, not is_white_turn))
    return None

solution = find_solution()
print(solution if solution else "No")