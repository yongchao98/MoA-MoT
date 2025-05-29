from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_player_pos(board):
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['*', '%']:
                return (i, j)
    return None

def is_goal(board):
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] == '@' or (board[i][j] in ['*', '%'] and board[i][j] != 'X'):
                return False
    return True

def move(board, direction):
    new_board = copy.deepcopy(board)
    player = get_player_pos(new_board)
    if not player:
        return None
    
    px, py = player
    dx, dy = {
        'U': (-1, 0),
        'D': (1, 0),
        'L': (0, -1),
        'R': (0, 1)
    }[direction]
    
    nx, ny = px + dx, py + dy
    if not is_valid(nx, ny, len(board), len(board[0])) or new_board[nx][ny] == '+':
        return None
    
    # Moving onto empty space or goal
    if new_board[nx][ny] in ['-', 'X']:
        new_board[px][py] = 'X' if board[px][py] == '%' else '-'
        new_board[nx][ny] = '%' if new_board[nx][ny] == 'X' else '*'
        return new_board
    
    # Pushing a box
    if new_board[nx][ny] in ['@', '$']:
        nnx, nny = nx + dx, ny + dy
        if not is_valid(nnx, nny, len(board), len(board[0])) or \
           new_board[nnx][nny] in ['+', '@', '$']:
            return None
        
        new_board[px][py] = 'X' if board[px][py] == '%' else '-'
        new_board[nx][ny] = '%' if board[nx][ny] == '$' else '*'
        new_board[nnx][nny] = '$' if new_board[nnx][nny] == 'X' else '@'
        return new_board
    
    return None

def board_to_string(board):
    return '\n'.join(''.join(row) for row in board)

def solve_sokoban(initial_board):
    queue = deque([(initial_board, "")])
    seen = {board_to_string(initial_board)}
    
    while queue:
        board, path = queue.popleft()
        
        if is_goal(board):
            return path
        
        for direction in ['U', 'D', 'L', 'R']:
            new_board = move(board, direction)
            if new_board:
                board_str = board_to_string(new_board)
                if board_str not in seen:
                    seen.add(board_str)
                    queue.append((new_board, path + direction))
    
    return None

# Initial board
initial_board = [
    ['+', '+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '-', '-', '-', '-', '-', '-', 'X', '+'],
    ['+', '@', '-', '-', '-', 'X', '$', '-', '+'],
    ['+', '%', '@', '-', '-', '-', '@', '-', '+'],
    ['+', '@', '-', '-', '-', '@', 'X', '-', '+'],
    ['+', '-', '@', '-', '-', '-', 'X', 'X', '+'],
    ['+', 'X', '+', '-', '-', '-', '-', '@', '+'],
    ['+', '-', '-', '-', '-', '-', '-', '-', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(initial_board)
print(f"<<<{solution}>>>")