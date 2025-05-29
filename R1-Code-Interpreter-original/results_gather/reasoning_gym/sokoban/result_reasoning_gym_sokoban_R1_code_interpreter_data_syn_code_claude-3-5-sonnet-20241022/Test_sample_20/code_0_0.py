from collections import deque
import copy

def is_valid(board, row, col):
    return 0 <= row < len(board) and 0 <= col < len(board[0]) and board[row][col] != '+'

def get_next_states(board, player_pos):
    moves = [('U', -1, 0), ('D', 1, 0), ('L', 0, -1), ('R', 0, 1)]
    next_states = []
    
    for move, dr, dc in moves:
        r, c = player_pos
        new_r, new_c = r + dr, c + dc
        
        if not is_valid(board, new_r, new_c):
            continue
            
        new_board = [list(row) for row in board]
        
        # If moving to empty space or goal
        if board[new_r][new_c] in '-X':
            new_board[r][c] = '-' if board[r][c] == '*' else 'X'
            new_board[new_r][new_c] = '*' if board[new_r][new_c] == '-' else '%'
            next_states.append((move, new_board, (new_r, new_c)))
            
        # If moving box
        elif board[new_r][new_c] in '@$':
            box_r, box_c = new_r + dr, new_c + dc
            if is_valid(board, box_r, box_c) and board[box_r][box_c] in '-X':
                new_board[r][c] = '-' if board[r][c] == '*' else 'X'
                new_board[new_r][new_c] = '*' if board[new_r][new_c] == '@' else '%'
                new_board[box_r][box_c] = '@' if board[box_r][box_c] == '-' else '$'
                next_states.append((move, new_board, (new_r, new_c)))
    
    return next_states

def board_to_string(board):
    return '\n'.join(''.join(row) for row in board)

def is_goal(board):
    for row in board:
        for cell in row:
            if cell == '@' or cell == 'X':
                return False
    return True

def solve_sokoban(initial_board):
    initial_state = [list(row) for row in initial_board.split('\n')]
    
    # Find player position
    player_pos = None
    for i in range(len(initial_state)):
        for j in range(len(initial_state[0])):
            if initial_state[i][j] in '*%':
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    visited = set()
    queue = deque([([], initial_state, player_pos)])
    
    while queue:
        path, current_board, pos = queue.popleft()
        board_str = board_to_string(current_board)
        
        if board_str in visited:
            continue
            
        visited.add(board_str)
        
        if is_goal(current_board):
            return ''.join(path)
            
        for move, new_board, new_pos in get_next_states(current_board, pos):
            queue.append((path + [move], new_board, new_pos))
    
    return None

# Initial board
initial_board = """++++++++
+XX@-@-+
+--XX@-+
+--$@-+
+-*---+
++++++++"""

solution = solve_sokoban(initial_board)
print(solution)