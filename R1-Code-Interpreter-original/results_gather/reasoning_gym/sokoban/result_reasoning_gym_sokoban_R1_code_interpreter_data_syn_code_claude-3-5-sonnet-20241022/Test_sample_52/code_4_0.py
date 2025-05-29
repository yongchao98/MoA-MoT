def parse_board(board_str):
    return [list(row.strip()) for row in board_str.strip().split('\n')]

def get_player_pos(board):
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['@', '*']:
                return (i, j)
    return None

def is_within_bounds(pos, board):
    return 0 <= pos[0] < len(board) and 0 <= pos[1] < len(board[0])

def is_valid_move(pos, direction, board):
    new_pos = (pos[0] + direction[0], pos[1] + direction[1])
    
    # Check boundaries and walls
    if not is_within_bounds(new_pos, board):
        return False
    if board[new_pos[0]][new_pos[1]] == '+':
        return False
        
    # If moving to a box position, check if box can be pushed
    if board[new_pos[0]][new_pos[1]] in ['$', '@']:
        box_new_pos = (new_pos[0] + direction[0], new_pos[1] + direction[1])
        if not is_within_bounds(box_new_pos, board):
            return False
        if board[box_new_pos[0]][box_new_pos[1]] in ['+', '$', '@']:
            return False
    return True

def make_move(board, pos, direction):
    new_board = [row[:] for row in board]
    new_pos = (pos[0] + direction[0], pos[1] + direction[1])
    
    # Update old player position
    if board[pos[0]][pos[1]] == '*':
        new_board[pos[0]][pos[1]] = 'X'
    else:
        new_board[pos[0]][pos[1]] = '-'
    
    # Handle box pushing
    if board[new_pos[0]][new_pos[1]] in ['$', '@']:
        box_new_pos = (new_pos[0] + direction[0], new_pos[1] + direction[1])
        if board[box_new_pos[0]][box_new_pos[1]] == 'X':
            new_board[box_new_pos[0]][box_new_pos[1]] = '$'
        else:
            new_board[box_new_pos[0]][box_new_pos[1]] = '@'
    
    # Update new player position
    if board[new_pos[0]][new_pos[1]] == 'X':
        new_board[new_pos[0]][new_pos[1]] = '*'
    else:
        new_board[new_pos[0]][new_pos[1]] = '@'
    
    return new_board

def is_solved(board):
    for row in board:
        for cell in row:
            if cell == 'X':  # If there's an unfilled goal
                return False
    return True

def board_to_string(board):
    return '\n'.join(''.join(row) for row in board)

def solve_sokoban(board):
    directions = {
        'U': (-1, 0),
        'D': (1, 0),
        'L': (0, -1),
        'R': (0, 1)
    }
    
    from collections import deque
    queue = deque([(board, '', get_player_pos(board))])
    visited = set()
    max_moves = 20
    
    while queue:
        current_board, path, player_pos = queue.popleft()
        
        if len(path) > max_moves:
            continue
            
        board_state = board_to_string(current_board)
        if board_state in visited:
            continue
        visited.add(board_state)
        
        if is_solved(current_board):
            return path
        
        for move, direction in directions.items():
            if is_valid_move(player_pos, direction, current_board):
                new_board = make_move(current_board, player_pos, direction)
                new_pos = (player_pos[0] + direction[0], player_pos[1] + direction[1])
                queue.append((new_board, path + move, new_pos))
    
    return None

# Initial board
board_str = """
+ + + + + + + +
+ - - - - - - +
+ - $ @ - $ - +
+ - X - @ $ - +
+ - - - X * - +
+ + + + + + + +
"""

initial_board = parse_board(board_str)
solution = solve_sokoban(initial_board)
print(solution)