from collections import deque

def is_valid(x, y, board):
    return 0 <= x < len(board) and 0 <= y < len(board[0]) and board[x][y] != '+'

def get_state_key(board):
    return ''.join(''.join(row) for row in board)

def can_push_box(board, player_x, player_y, box_x, box_y, dx, dy):
    # Check if box can be pushed to new position
    new_box_x, new_box_y = box_x + dx, box_y + dy
    return (is_valid(new_box_x, new_box_y, board) and 
            board[new_box_x][new_box_y] in ['-', 'X'])

def solve_sokoban(initial):
    board = [list(row) for row in initial]
    rows, cols = len(board), len(board[0])
    
    # Find initial positions
    players = []
    boxes = []
    goals = []
    
    for i in range(rows):
        for j in range(cols):
            if board[i][j] == '@':
                players.append((i, j))
            elif board[i][j] == '*':
                players.append((i, j))
                goals.append((i, j))
            elif board[i][j] == '$':
                boxes.append((i, j))
            elif board[i][j] == 'X':
                goals.append((i, j))
    
    directions = {
        'U': (-1, 0),
        'D': (1, 0),
        'L': (0, -1),
        'R': (0, 1)
    }
    
    # Start with first player and first box
    start_player = players[1]  # Using the second player as starting point
    start_box = boxes[0]
    
    visited = set()
    queue = deque([(board, "", start_player, start_box)])
    
    while queue:
        current_board, moves, (px, py), (bx, by) = queue.popleft()
        
        # Check if box is on any goal
        if (bx, by) in goals:
            return moves
        
        state = get_state_key(current_board)
        if state in visited:
            continue
        visited.add(state)
        
        for move, (dx, dy) in directions.items():
            new_px, new_py = px + dx, py + dy
            
            if not is_valid(new_px, new_py, current_board):
                continue
                
            new_board = [row[:] for row in current_board]
            new_box_pos = (bx, by)
            
            # If moving to box position
            if (new_px, new_py) == (bx, by):
                if can_push_box(current_board, px, py, bx, by, dx, dy):
                    new_bx, new_by = bx + dx, by + dy
                    new_box_pos = (new_bx, new_by)
                    
                    # Update box position
                    new_board[bx][by] = '@'  # Player moves to box position
                    if (new_bx, new_by) in goals:
                        new_board[new_bx][new_by] = '$'
                    else:
                        new_board[new_bx][new_by] = '@'
                else:
                    continue
            else:
                # Regular player movement
                if current_board[new_px][new_py] in ['-', 'X']:
                    new_board[new_px][new_py] = '@'
            
            new_board[px][py] = '-'
            queue.append((new_board, moves + move, (new_px, new_py), new_box_pos))
    
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