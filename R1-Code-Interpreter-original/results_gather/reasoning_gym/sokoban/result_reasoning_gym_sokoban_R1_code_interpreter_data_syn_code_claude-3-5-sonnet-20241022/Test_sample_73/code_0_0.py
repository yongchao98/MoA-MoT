from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_player_pos(board):
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['*', '@', '%']:
                return (i, j)
    return None

def is_goal(board):
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] == '@' or (board[i][j] == 'X' and not any(b == '$' for row in board for b in row)):
                return False
    return True

def get_board_state(board):
    return '\n'.join(''.join(row) for row in board)

def solve_sokoban(initial_board):
    directions = [(0, 1, 'R'), (0, -1, 'L'), (1, 0, 'D'), (-1, 0, 'U')]
    queue = deque([(initial_board, "")])
    visited = set()
    
    while queue:
        current_board, path = queue.popleft()
        
        if is_goal(current_board):
            return path
            
        board_state = get_board_state(current_board)
        if board_state in visited:
            continue
        visited.add(board_state)
        
        player_pos = get_player_pos(current_board)
        if not player_pos:
            continue
            
        px, py = player_pos
        
        for dx, dy, move in directions:
            new_px, new_py = px + dx, py + dy
            
            if not is_valid(new_px, new_py, len(current_board), len(current_board[0])):
                continue
                
            if current_board[new_px][new_py] in ['+']:
                continue
                
            new_board = [list(row) for row in current_board]
            
            # If next position is a box
            if current_board[new_px][new_py] in ['@', '$']:
                box_px, box_py = new_px + dx, new_py + dy
                
                if not is_valid(box_px, box_py, len(current_board), len(current_board[0])):
                    continue
                    
                if current_board[box_px][box_py] in ['+', '@', '$']:
                    continue
                
                # Move box
                if current_board[box_px][box_py] == 'X':
                    new_board[box_px][box_py] = '$'
                else:
                    new_board[box_px][box_py] = '@'
                    
                if current_board[new_px][new_py] == '$':
                    new_board[new_px][new_py] = 'X'
                else:
                    new_board[new_px][new_py] = '-'
            
            # Move player
            if current_board[px][py] == '%':
                new_board[px][py] = 'X'
            else:
                new_board[px][py] = '-'
                
            if current_board[new_px][new_py] == 'X':
                new_board[new_px][new_py] = '%'
            else:
                new_board[new_px][new_py] = '*'
            
            queue.append((new_board, path + move))
    
    return None

# Initialize the board
initial_board = [
    ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '-', '-', '-', '-', '-', '-', '-', '-', '+'],
    ['+', '-', '$', '@', '-', '-', '@', '+', '+', '+'],
    ['+', '+', '+', 'X', '+', '-', 'X', '+', '+', '+'],
    ['+', '+', '+', '-', '-', '-', '-', '$', 'X', '+'],
    ['+', '+', '+', '-', '-', '-', '-', '@', '@', '+'],
    ['+', '+', '+', '$', 'X', '+', '+', '$', '*', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(initial_board)
print(solution if solution else "No solution found")