def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def get_state_key(board):
    boxes = []
    player = None
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['@', '$']:
                boxes.append((i, j))
            elif board[i][j] in ['*', '%']:
                player = (i, j)
    boxes.sort()  # Sort for consistent state representation
    return str(player) + '|' + str(boxes)

def is_deadlock(board, box_pos):
    rows, cols = len(board), len(board[0])
    x, y = box_pos
    
    # Corner deadlock
    if (board[x-1][y] == '+' and board[x][y-1] == '+') or \
       (board[x-1][y] == '+' and board[x][y+1] == '+') or \
       (board[x+1][y] == '+' and board[x][y-1] == '+') or \
       (board[x+1][y] == '+' and board[x][y+1] == '+'):
        return True
    return False

def solve_sokoban(board):
    rows, cols = len(board), len(board[0])
    goals = []
    boxes = []
    player = None
    
    # Find initial positions
    for i in range(rows):
        for j in range(cols):
            if board[i][j] in ['X', '%', '$']:
                goals.append((i, j))
            if board[i][j] in ['@', '$']:
                boxes.append((i, j))
            if board[i][j] in ['*', '%']:
                player = (i, j)
    
    def get_moves(pos, board):
        moves = []
        for dx, dy, move in [(-1,0,'U'), (1,0,'D'), (0,-1,'L'), (0,1,'R')]:
            new_x, new_y = pos[0] + dx, pos[1] + dy
            
            if board[new_x][new_y] == '+':
                continue
                
            if board[new_x][new_y] in ['@', '$']:
                box_x, box_y = new_x + dx, new_y + dy
                if board[box_x][box_y] not in ['+', '@', '$'] and not is_deadlock(board, (box_x, box_y)):
                    moves.append((move, (new_x, new_y), (box_x, box_y)))
            else:
                moves.append((move, (new_x, new_y), None))
        return moves
    
    from collections import deque
    queue = deque([(board, player, "", 0)])
    visited = set()
    max_depth = 50  # Limit search depth
    
    while queue:
        current_board, pos, path, depth = queue.popleft()
        
        if depth > max_depth:
            continue
            
        state = get_state_key(current_board)
        if state in visited:
            continue
        visited.add(state)
        
        # Check win condition
        all_boxes_on_goals = True
        for i in range(rows):
            for j in range(cols):
                if current_board[i][j] == '@':
                    all_boxes_on_goals = False
                    break
            if not all_boxes_on_goals:
                break
        
        if all_boxes_on_goals:
            return path
        
        # Get valid moves
        for move, new_pos, box_pos in get_moves(pos, current_board):
            new_board = [list(row) for row in current_board]
            
            # Update player
            new_board[pos[0]][pos[1]] = '-' if current_board[pos[0]][pos[1]] == '*' else 'X'
            new_board[new_pos[0]][new_pos[1]] = '*'
            
            # Update box if pushing
            if box_pos:
                new_board[new_pos[0]][new_pos[1]] = '-'
                new_board[box_pos[0]][box_pos[1]] = '@'
                if (box_pos[0], box_pos[1]) in goals:
                    new_board[box_pos[0]][box_pos[1]] = '$'
            
            queue.append((new_board, new_pos, path + move, depth + 1))
    
    return None

# Initialize board
board = [
    ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '-', '@', '-', 'X', '-', '-', '-', '-', '+'],
    ['+', '-', 'X', '-', '-', '@', 'X', '@', '-', '+'],
    ['+', '-', '-', '*', '$', '-', '-', '-', 'X', '+'],
    ['+', '-', '$', '-', '@', '-', '-', '-', '-', '+'],
    ['+', '-', '-', '-', 'X', '@', '-', '-', '-', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(board)
print(solution if solution else "No solution found")