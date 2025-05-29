from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_pos(x, y, dx, dy):
    return x + dx, y + dy

def get_state_key(board):
    return ''.join(''.join(row) for row in board)

def is_win(board):
    for row in board:
        for cell in row:
            if cell == '@' or cell == 'X':  # If there's a box not on goal or goal without box
                return False
    return True

def get_moves(board, player_pos):
    rows, cols = len(board), len(board[0])
    moves = []
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    
    for move, (dx, dy) in directions.items():
        new_x, new_y = get_next_pos(player_pos[0], player_pos[1], dx, dy)
        
        if not is_valid(new_x, new_y, rows, cols) or board[new_x][new_y] == '+':
            continue
            
        if board[new_x][new_y] in ['@', '$']:
            box_x, box_y = get_next_pos(new_x, new_y, dx, dy)
            if not is_valid(box_x, box_y, rows, cols) or board[box_x][box_y] in ['+', '@', '$']:
                continue
            moves.append((move, (new_x, new_y), (box_x, box_y)))
        else:
            moves.append((move, (new_x, new_y), None))
    
    return moves

def solve_sokoban(initial_board):
    rows, cols = len(initial_board), len(initial_board[0])
    player_pos = None
    
    # Find player position
    for i in range(rows):
        for j in range(cols):
            if initial_board[i][j] in ['*', '%']:
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    queue = deque([(initial_board, player_pos, "")])
    visited = set()
    
    while queue:
        board, pos, path = queue.popleft()
        state_key = get_state_key(board)
        
        if state_key in visited:
            continue
            
        visited.add(state_key)
        
        if is_win(board):
            return path
            
        for move, new_pos, box_pos in get_moves(board, pos):
            new_board = [list(row) for row in board]
            
            # Update player position
            old_cell = 'X' if board[pos[0]][pos[1]] == '%' else '-'
            new_board[pos[0]][pos[1]] = old_cell
            
            new_cell = '%' if board[new_pos[0]][new_pos[1]] == 'X' else '*'
            new_board[new_pos[0]][new_pos[1]] = new_cell
            
            # Update box position if pushing
            if box_pos:
                old_box_cell = board[new_pos[0]][new_pos[1]]
                new_box_cell = '$' if board[box_pos[0]][box_pos[1]] == 'X' else '@'
                new_board[box_pos[0]][box_pos[1]] = new_box_cell
            
            queue.append((new_board, new_pos, path + move))
    
    return None

# Initialize the board
initial_board = [
    ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '-', '@', '-', 'X', '-', '-', '-', '-', '+'],
    ['+', '-', 'X', '-', '-', '@', 'X', '@', '-', '+'],
    ['+', '-', '-', '*', '$', '-', '-', '-', 'X', '+'],
    ['+', '-', '$', '-', '@', '-', '-', '-', '-', '+'],
    ['+', '-', '-', '-', 'X', '@', '-', '-', '-', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(initial_board)
print(solution if solution else "No solution found")