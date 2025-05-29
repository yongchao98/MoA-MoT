from collections import deque

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_state_key(board, player_pos):
    state = ''
    for row in board:
        state += ''.join(row)
    return state + f',{player_pos[0]},{player_pos[1]}'

def is_win(board):
    goals = sum(row.count('X') for row in board)
    boxes_on_goals = sum(row.count('$') for row in board)
    return goals == 0 and boxes_on_goals == 3

def get_moves(board, player_pos):
    rows, cols = len(board), len(board[0])
    moves = []
    directions = [(-1,0,'U'), (1,0,'D'), (0,-1,'L'), (0,1,'R')]
    
    for dx, dy, move in directions:
        new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
        
        if not is_valid(new_x, new_y, rows, cols) or board[new_x][new_y] == '+':
            continue
            
        if board[new_x][new_y] in ['-', 'X']:
            moves.append(((new_x, new_y), move, False))  # Fixed tuple structure
        elif board[new_x][new_y] in ['@', '$']:
            push_x, push_y = new_x + dx, new_y + dy
            if (is_valid(push_x, push_y, rows, cols) and 
                board[push_x][push_y] in ['-', 'X']):
                moves.append(((new_x, new_y), move, True))  # Fixed tuple structure
    
    return moves

def make_move(board, player_pos, new_pos, is_push):
    new_board = [row[:] for row in board]
    px, py = player_pos
    nx, ny = new_pos
    
    # Clear old player position
    if board[px][py] == '%':
        new_board[px][py] = 'X'
    else:
        new_board[px][py] = '-'
    
    # Handle box pushing
    if is_push:
        dx, dy = nx - px, ny - py
        push_x, push_y = nx + dx, ny + dy
        
        if new_board[nx][ny] == '@':
            if new_board[push_x][push_y] == 'X':
                new_board[push_x][push_y] = '$'
            else:
                new_board[push_x][push_y] = '@'
            new_board[nx][ny] = '-'
        else:  # box on goal
            if new_board[push_x][push_y] == 'X':
                new_board[push_x][push_y] = '$'
            else:
                new_board[push_x][push_y] = '@'
            new_board[nx][ny] = 'X'
    
    # Place player in new position
    if new_board[nx][ny] == 'X':
        new_board[nx][ny] = '%'
    else:
        new_board[nx][ny] = '*'
    
    return new_board

def solve_sokoban(board):
    rows, cols = len(board), len(board[0])
    player_pos = None
    
    # Find initial player position
    for i in range(rows):
        for j in range(cols):
            if board[i][j] in ['*', '%']:
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    queue = deque([(board, player_pos, "")])
    visited = {get_state_key(board, player_pos)}
    
    while queue:
        current_board, current_pos, path = queue.popleft()
        
        if is_win(current_board):
            return path
        
        for new_pos, move, is_push in get_moves(current_board, current_pos):
            new_board = make_move(current_board, current_pos, new_pos, is_push)
            state_key = get_state_key(new_board, new_pos)
            
            if state_key not in visited:
                visited.add(state_key)
                queue.append((new_board, new_pos, path + move))
    
    return None

# Initialize the puzzle
initial_board = [
    list('++++++'),
    list('+XX@-+'),
    list('+-X--+'),
    list('+@-@-+'),
    list('+%@--+'),
    list('++++++')
]

solution = solve_sokoban(initial_board)
print(solution)