from collections import deque
import copy

def is_valid(x, y, rows, cols, board):
    return 0 <= x < rows and 0 <= y < cols and board[x][y] != '+'

def get_player_pos(board):
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] == '%':
                return (i, j)
    return None

def get_state_key(board):
    # Create a unique key for the board state
    state = []
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] in ['@', '$', '%']:
                state.append((i, j, board[i][j]))
    return tuple(sorted(state))

def solve_sokoban(initial_board):
    moves = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    queue = deque([(initial_board, "")])
    visited = set()
    max_depth = 50  # Limit search depth
    
    while queue and len(queue[0][1]) < max_depth:
        board, path = queue.popleft()
        state_key = get_state_key(board)
        
        if state_key in visited:
            continue
        visited.add(state_key)
        
        # Check if all boxes are on goals
        boxes_on_goals = sum(1 for row in board for cell in row if cell == '$')
        total_boxes = sum(1 for row in board for cell in row if cell in ['@', '$'])
        if boxes_on_goals == total_boxes:
            return path
        
        px, py = get_player_pos(board)
        if px is None:
            continue
        
        for move, (dx, dy) in moves.items():
            new_px, new_py = px + dx, py + dy
            
            if not is_valid(new_px, new_py, len(board), len(board[0]), board):
                continue
            
            new_board = [list(row) for row in board]
            
            # Handle box pushing
            if new_board[new_px][new_py] in ['@', '$']:
                box_px, box_py = new_px + dx, new_py + dy
                
                if not is_valid(box_px, box_py, len(board), len(board[0]), board):
                    continue
                if new_board[box_px][box_py] in ['@', '$', '+']:
                    continue
                
                # Move box
                is_box_on_goal = new_board[box_px][box_py] == 'X'
                new_board[box_px][box_py] = '$' if is_box_on_goal else '@'
                new_board[new_px][new_py] = '%'
            else:
                # Move player
                is_player_on_goal = new_board[new_px][new_py] == 'X'
                new_board[new_px][new_py] = '%' if is_player_on_goal else '%'
            
            # Update old player position
            new_board[px][py] = 'X' if board[px][py] == '%' else '-'
            
            queue.append((new_board, path + move))
    
    return None

# Initial board
initial_board = [
    ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '-', '-', '-', '-', '-', '-', '-', '+', '+'],
    ['+', '-', '-', '@', '-', '-', '-', '-', '+', '+'],
    ['+', '-', '-', '-', 'X', 'X', '@', '-', '-', '+'],
    ['+', '-', '$', '$', '$', '-', '@', '-', '@', '+'],
    ['+', '-', '-', '-', '-', '-', '-', 'X', '%', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(initial_board)
print(solution)