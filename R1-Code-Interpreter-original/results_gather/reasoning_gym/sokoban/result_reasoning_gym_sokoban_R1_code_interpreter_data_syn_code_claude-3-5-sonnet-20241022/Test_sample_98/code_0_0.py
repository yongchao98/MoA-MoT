from collections import deque
import copy

def is_valid_pos(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_pos(x, y, direction):
    if direction == 'U': return (x-1, y)
    if direction == 'D': return (x+1, y)
    if direction == 'L': return (x, y-1)
    if direction == 'R': return (x, y+1)

def get_state_key(board, player_pos):
    # Convert current state to a string for hashing
    return (''.join([''.join(row) for row in board]), player_pos)

def is_win(board):
    # Check if all boxes are on goals
    box_count = sum(row.count('@') for row in board)
    box_on_goal = sum(row.count('$') for row in board)
    return box_count == 0 and box_on_goal > 0

def solve_sokoban(initial_board):
    rows, cols = len(initial_board), len(initial_board[0])
    board = [list(row) for row in initial_board]
    
    # Find player position
    player_pos = None
    for i in range(rows):
        for j in range(cols):
            if board[i][j] == '*':
                player_pos = (i, j)
                if any(board[x][y] == 'X' for x, y in [(i,j)]):
                    board[i][j] = '%'
                else:
                    board[i][j] = '-'
                break
        if player_pos: break

    visited = set()
    queue = deque([(board, player_pos, "")])
    
    while queue:
        current_board, current_pos, path = queue.popleft()
        state_key = get_state_key(current_board, current_pos)
        
        if state_key in visited:
            continue
        
        visited.add(state_key)
        
        if is_win(current_board):
            return path
        
        for direction in ['U', 'D', 'L', 'R']:
            new_pos = get_next_pos(current_pos[0], current_pos[1], direction)
            
            if not is_valid_pos(new_pos[0], new_pos[1], rows, cols):
                continue
                
            new_board = [row[:] for row in current_board]
            
            # Check if move is valid
            if new_board[new_pos[0]][new_pos[1]] in ['+']:
                continue
                
            if new_board[new_pos[0]][new_pos[1]] in ['@', '$']:
                # There's a box, check if we can push it
                box_new_pos = get_next_pos(new_pos[0], new_pos[1], direction)
                
                if not is_valid_pos(box_new_pos[0], box_new_pos[1], rows, cols):
                    continue
                    
                if new_board[box_new_pos[0]][box_new_pos[1]] in ['+', '@', '$']:
                    continue
                
                # Move the box
                if new_board[box_new_pos[0]][box_new_pos[1]] == 'X':
                    new_board[box_new_pos[0]][box_new_pos[1]] = '$'
                else:
                    new_board[box_new_pos[0]][box_new_pos[1]] = '@'
                
                new_board[new_pos[0]][new_pos[1]] = '-'
            
            queue.append((new_board, new_pos, path + direction))
    
    return None

# Initial board
board = [
    "++++++++",
    "+"----X+",
    "+-@X@--+",
    "+X-*X--+",
    "+-@XX@-+",
    "+@--@--+",
    "+------+",
    "++++++++"
]

solution = solve_sokoban(board)
print(solution)