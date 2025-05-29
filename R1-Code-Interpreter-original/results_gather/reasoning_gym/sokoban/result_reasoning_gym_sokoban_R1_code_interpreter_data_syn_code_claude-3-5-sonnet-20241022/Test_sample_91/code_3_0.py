from collections import deque

def is_valid(x, y, board):
    return 0 <= x < len(board) and 0 <= y < len(board[0]) and board[x][y] != '+'

def get_state_string(board):
    return ''.join(''.join(row) for row in board)

def solve_sokoban(initial):
    board = [list(row) for row in initial]
    rows, cols = len(board), len(board[0])
    
    # Find active player (starting with topmost, then leftmost)
    player_pos = None
    for i in range(rows):
        for j in range(cols):
            if board[i][j] in ['@', '*']:
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    directions = [(-1, 0, 'U'), (1, 0, 'D'), (0, -1, 'L'), (0, 1, 'R')]
    visited = set()
    queue = deque([(board, "", player_pos)])
    
    while queue:
        current_board, moves, (px, py) = queue.popleft()
        
        # Check if won (all boxes on goals)
        won = True
        for i in range(rows):
            for j in range(cols):
                if current_board[i][j] == '@':  # Box not on goal
                    won = False
                    break
            if not won:
                break
        if won:
            return moves
        
        state = get_state_string(current_board)
        if state in visited:
            continue
        visited.add(state)
        
        # Try each direction
        for dx, dy, move in directions:
            new_x, new_y = px + dx, py + dy
            
            if not is_valid(new_x, new_y, current_board):
                continue
            
            new_board = [row[:] for row in current_board]
            can_move = False
            
            # If moving to empty space or goal
            if current_board[new_x][new_y] in ['-', 'X']:
                can_move = True
                # Update player position
                if current_board[px][py] == '*':
                    new_board[px][py] = 'X'
                else:
                    new_board[px][py] = '-'
                
                if current_board[new_x][new_y] == 'X':
                    new_board[new_x][new_y] = '*'
                else:
                    new_board[new_x][new_y] = '@'
            
            # If pushing a box
            elif current_board[new_x][new_y] == '@':
                push_x, push_y = new_x + dx, new_y + dy
                
                if (is_valid(push_x, push_y, current_board) and 
                    current_board[push_x][push_y] in ['-', 'X']):
                    can_move = True
                    # Move box
                    if current_board[push_x][push_y] == 'X':
                        new_board[push_x][push_y] = '$'
                    else:
                        new_board[push_x][push_y] = '@'
                    
                    # Move player
                    new_board[new_x][new_y] = '@'
                    new_board[px][py] = '-'
            
            if can_move:
                queue.append((new_board, moves + move, (new_x, new_y)))
    
    return None

initial = [
    "++++++",
    "+-@-X+",
    "+-$-++",
    "+X@X-+",
    "+-@--+",
    "+-@X*+",
    "+----+",
    "++++++"
]

solution = solve_sokoban(initial)
print(solution)