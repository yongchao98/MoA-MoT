from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_states(board, player_pos):
    rows, cols = len(board), len(board[0])
    moves = [(-1, 0, 'U'), (1, 0, 'D'), (0, -1, 'L'), (0, 1, 'R')]
    next_states = []
    
    for dx, dy, move in moves:
        new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
        
        if not is_valid(new_x, new_y, rows, cols) or board[new_x][new_y] == '+':
            continue
            
        new_board = [row[:] for row in board]
        new_player_pos = [new_x, new_y]
        
        if board[new_x][new_y] in ['@', '$']:
            box_x, box_y = new_x + dx, new_y + dy
            if not is_valid(box_x, box_y, rows, cols) or board[box_x][box_y] in ['+', '@', '$']:
                continue
                
            # Move box
            if board[new_x][new_y] == '@':
                new_board[new_x][new_y] = '-'
                new_board[box_x][box_y] = '@' if board[box_x][box_y] == '-' else '$'
            else:  # box on goal
                new_board[new_x][new_y] = 'X'
                new_board[box_x][box_y] = '@' if board[box_x][box_y] == '-' else '$'
        
        # Move player
        if board[player_pos[0]][player_pos[1]] == '*':
            new_board[player_pos[0]][player_pos[1]] = '-'
        else:  # player on goal
            new_board[player_pos[0]][player_pos[1]] = 'X'
            
        new_board[new_x][new_y] = '*' if board[new_x][new_y] in ['-', '@'] else '%'
        
        next_states.append((new_board, new_player_pos, move))
    return next_states

def is_solved(board):
    for row in board:
        for cell in row:
            if cell == '@' or cell == 'X':
                return False
    return True

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_sokoban(initial_board):
    rows, cols = len(initial_board), len(initial_board[0])
    player_pos = None
    
    # Find player position
    for i in range(rows):
        for j in range(cols):
            if initial_board[i][j] in ['*', '%']:
                player_pos = [i, j]
                break
        if player_pos:
            break
    
    queue = deque([(initial_board, player_pos, "")])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, current_pos, path = queue.popleft()
        
        if is_solved(current_board):
            return path
        
        for new_board, new_pos, move in get_next_states(current_board, current_pos):
            board_str = board_to_string(new_board)
            if board_str not in visited:
                visited.add(board_str)
                queue.append((new_board, new_pos, path + move))
    
    return None

# Initialize the board
initial_board = [
    list("+++++++++"),
    list("+*@---X$+"),
    list("+@-@-X+++"),
    list("+---$$$+"),
    list("+-++++++++"),
    list("+--+----+"),
    list("+X@--+++"),
    list("++-@X+++"),
    list("++X++++"),
    list("+++++++++")
]

solution = solve_sokoban(initial_board)
print(solution if solution else "No solution found")