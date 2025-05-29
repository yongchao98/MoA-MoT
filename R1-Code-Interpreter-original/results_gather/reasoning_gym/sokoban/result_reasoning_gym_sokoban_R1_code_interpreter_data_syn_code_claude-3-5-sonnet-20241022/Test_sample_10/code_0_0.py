from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_state_key(board):
    return ''.join([''.join(row) for row in board])

def is_win(board):
    box_count = sum(row.count('@') for row in board)
    goal_box_count = sum(row.count('$') for row in board)
    return box_count == 0 and goal_box_count == 3

def get_moves(board, player_pos):
    rows, cols = len(board), len(board[0])
    moves = []
    directions = [(-1,0,'U'), (1,0,'D'), (0,-1,'L'), (0,1,'R')]
    
    for dx, dy, move in directions:
        new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
        
        if not is_valid(new_x, new_y, rows, cols) or board[new_x][new_y] == '+':
            continue
            
        if board[new_x][new_y] in ['-', 'X']:
            moves.append((new_x, new_y, move, False))
        elif board[new_x][new_y] in ['@', '$']:
            push_x, push_y = new_x + dx, new_y + dy
            if (is_valid(push_x, push_y, rows, cols) and 
                board[push_x][push_y] in ['-', 'X']):
                moves.append((new_x, new_y, move, True))
                
    return moves

def solve_sokoban(initial_board):
    rows, cols = len(initial_board), len(initial_board[0])
    queue = deque()
    visited = set()
    
    # Find initial player position
    player_pos = None
    for i in range(rows):
        for j in range(cols):
            if initial_board[i][j] in ['*', '%']:
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    queue.append((initial_board, player_pos, ""))
    visited.add(get_state_key(initial_board))
    
    while queue:
        current_board, current_pos, path = queue.popleft()
        
        if is_win(current_board):
            return path
        
        moves = get_moves(current_board, current_pos)
        for new_x, new_y, move, is_push in moves:
            new_board = [list(row) for row in current_board]
            px, py = current_pos
            
            # Update player position
            if new_board[px][py] == '%':
                new_board[px][py] = 'X'
            else:
                new_board[px][py] = '-'
                
            # Handle box pushing
            if is_push:
                dx, dy = new_x - px, new_y - py
                push_x, push_y = new_x + dx, new_y + dy
                
                if new_board[new_x][new_y] == '@':
                    new_board[new_x][new_y] = '-'
                else:  # '$'
                    new_board[new_x][new_y] = 'X'
                    
                if new_board[push_x][push_y] == 'X':
                    new_board[push_x][push_y] = '$'
                else:
                    new_board[push_x][push_y] = '@'
            
            # Update player's new position
            if new_board[new_x][new_y] == 'X':
                new_board[new_x][new_y] = '%'
            else:
                new_board[new_x][new_y] = '*'
            
            state_key = get_state_key(new_board)
            if state_key not in visited:
                visited.add(state_key)
                queue.append((new_board, (new_x, new_y), path + move))
    
    return None

# Initialize the puzzle
initial_board = [
    ['+', '+', '+', '+', '+', '+'],
    ['+', 'X', 'X', '@', '-', '+'],
    ['+', '-', 'X', '-', '-', '+'],
    ['+', '@', '-', '@', '-', '+'],
    ['+', '%', '@', '-', '-', '+'],
    ['+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(initial_board)
print(solution)