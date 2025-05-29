from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_states(board, rows, cols):
    # Find player position
    player_x, player_y = None, None
    for i in range(rows):
        for j in range(cols):
            if board[i][j] in ['*', '%']:
                player_x, player_y = i, j
                break
        if player_x is not None:
            break

    moves = [('U', -1, 0), ('D', 1, 0), ('L', 0, -1), ('R', 0, 1)]
    next_states = []
    
    for move, dx, dy in moves:
        new_x, new_y = player_x + dx, player_y + dy
        
        if not is_valid(new_x, new_y, rows, cols) or board[new_x][new_y] == '+':
            continue
            
        new_board = [row[:] for row in board]
        
        # Check if moving onto/pushing a box
        if board[new_x][new_y] in ['@', '$']:
            box_x, box_y = new_x + dx, new_y + dy
            
            if not is_valid(box_x, box_y, rows, cols) or \
               board[box_x][box_y] in ['+', '@', '$']:
                continue
                
            # Move box
            new_board[box_x][box_y] = '$' if board[box_x][box_y] == 'X' else '@'
            new_board[new_x][new_y] = '%' if board[new_x][new_y] == '$' else '*'
        else:
            new_board[new_x][new_y] = '%' if board[new_x][new_y] == 'X' else '*'
            
        # Update old player position
        new_board[player_x][player_y] = 'X' if board[player_x][player_y] == '%' else '-'
        
        next_states.append((new_board, move))
    
    return next_states

def is_solved(board):
    return all(cell != 'X' for row in board for cell in row)

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_sokoban(initial_board):
    rows, cols = len(initial_board), len(initial_board[0])
    queue = deque([(initial_board, "")])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, path = queue.popleft()
        
        if is_solved(current_board):
            return path
            
        for next_board, move in get_next_states(current_board, rows, cols):
            board_str = board_to_string(next_board)
            if board_str not in visited:
                visited.add(board_str)
                queue.append((next_board, path + move))
    
    return None

# Initialize the board
initial_board = [
    list("+++++++" ),
    list("+++++X+" ),
    list("++++-@+" ),
    list("+++-@*+" ),
    list("++$-X$+" ),
    list("+-+-+++" ),
    list("+@+++++" ),
    list("+X+++++" ),
    list("+++$+$+" ),
    list("+++++++" )
]

solution = solve_sokoban(initial_board)
print(solution if solution else "No solution found")