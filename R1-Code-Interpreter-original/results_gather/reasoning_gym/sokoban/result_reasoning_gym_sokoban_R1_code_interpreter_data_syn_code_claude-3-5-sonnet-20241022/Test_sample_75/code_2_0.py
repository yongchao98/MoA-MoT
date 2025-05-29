from collections import deque

def get_valid_moves(board, pos):
    moves = []
    directions = {
        'U': (-1, 0),
        'D': (1, 0),
        'L': (0, -1),
        'R': (0, 1)
    }
    
    for move, (dx, dy) in directions.items():
        new_x, new_y = pos[0] + dx, pos[1] + dy
        
        # Check if move is within bounds
        if board[new_x][new_y] == '+':
            continue
            
        # If moving to empty space or goal
        if board[new_x][new_y] in ['-', 'X']:
            moves.append(move)
            
        # If moving to box, check if we can push it
        elif board[new_x][new_y] in ['$', '@']:
            box_x, box_y = new_x + dx, new_y + dy
            if board[box_x][box_y] in ['-', 'X']:
                moves.append(move)
                
    return moves

def make_move(board, pos, move):
    new_board = [list(row) for row in board]
    x, y = pos
    
    directions = {
        'U': (-1, 0),
        'D': (1, 0),
        'L': (0, -1),
        'R': (0, 1)
    }
    
    dx, dy = directions[move]
    new_x, new_y = x + dx, y + dy
    
    # Moving to empty space
    if new_board[new_x][new_y] in ['-', 'X']:
        # Update player position
        new_board[x][y] = '-'
        new_board[new_x][new_y] = '@'
        return [''.join(row) for row in new_board], (new_x, new_y)
        
    # Moving box
    elif new_board[new_x][new_y] in ['$', '@']:
        box_x, box_y = new_x + dx, new_y + dy
        # Update box position
        new_board[box_x][box_y] = '@' if new_board[box_x][box_y] == '-' else '$'
        # Update player position
        new_board[x][y] = '-'
        new_board[new_x][new_y] = '@'
        return [''.join(row) for row in new_board], (new_x, new_y)

def is_solved(board):
    goals = 0
    boxes_on_goals = 0
    for row in board:
        for cell in row:
            if cell == 'X':
                goals += 1
            elif cell == '$':
                boxes_on_goals += 1
    return goals == 0 and boxes_on_goals == 3  # We have 3 goals to fill

def solve_sokoban(board):
    # Find initial player position (using the bottom-most player)
    player_pos = None
    for i in range(len(board)-1, -1, -1):
        for j in range(len(board[i])):
            if board[i][j] in ['@', '*']:
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    # Clean up board - remove other player positions
    board = [list(row) for row in board]
    for i in range(len(board)):
        for j in range(len(board[i])):
            if (i, j) != player_pos and board[i][j] in ['@', '*']:
                board[i][j] = '-'
    board = [''.join(row) for row in board]
    
    # BFS
    queue = deque([(board, player_pos, "")])
    seen = {'\n'.join(board)}
    
    while queue:
        current_board, pos, path = queue.popleft()
        
        if is_solved(current_board):
            return path
            
        for move in get_valid_moves(current_board, pos):
            new_board, new_pos = make_move(current_board, pos, move)
            board_key = '\n'.join(new_board)
            
            if board_key not in seen:
                seen.add(board_key)
                queue.append((new_board, new_pos, path + move))
    
    return None

# Initial state
initial_state = [
    "+++++++" ,
    "+++$+++",
    "+++$+++",
    "++++$$+",
    "+-----+",
    "+----X+",
    "+-X-@-+",
    "+---@*+",
    "+--X@-+",
    "+++++++"
]

solution = solve_sokoban(initial_state)
print(solution)