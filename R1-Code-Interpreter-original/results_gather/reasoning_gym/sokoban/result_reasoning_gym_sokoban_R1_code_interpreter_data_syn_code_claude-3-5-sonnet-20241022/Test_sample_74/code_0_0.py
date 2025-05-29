from collections import deque

def get_width(board):
    return board.index('\n')

def get_pos(board, char):
    pos = []
    for i, c in enumerate(board):
        if c in char:
            pos.append((i // (get_width(board) + 1), i % (get_width(board) + 1)))
    return pos

def is_valid_move(board, player_pos, dx, dy):
    width = get_width(board)
    new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
    new_pos = new_x * (width + 1) + new_y
    
    if board[new_pos] == '+':
        return False, None
    
    if board[new_pos] in '@$':
        box_x, box_y = new_x + dx, new_y + dy
        box_pos = box_x * (width + 1) + box_y
        if board[box_pos] in '+@$':
            return False, None
        return True, box_pos
    
    return True, None

def make_move(board, player_pos, dx, dy):
    width = get_width(board)
    pos = player_pos[0] * (width + 1) + player_pos[1]
    new_x, new_y = player_pos[0] + dx, player_pos[1] + dy
    new_pos = new_x * (width + 1) + new_y
    
    board = list(board)
    
    # Update player position
    if board[pos] == '%':
        board[pos] = 'X'
    else:
        board[pos] = '-'
        
    if board[new_pos] == 'X':
        board[new_pos] = '%'
    else:
        board[new_pos] = '*'
    
    # If pushing a box
    if board[new_pos] in '@$':
        box_x, box_y = new_x + dx, new_y + dy
        box_pos = box_x * (width + 1) + box_y
        
        if board[box_pos] == 'X':
            board[box_pos] = '$'
        else:
            board[box_pos] = '@'
            
        if board[new_pos] == '$':
            board[new_pos] = '%'
        else:
            board[new_pos] = '*'
    
    return ''.join(board)

def solve_sokoban(initial_board):
    moves = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    queue = deque([(initial_board, "")])
    visited = {initial_board}
    
    while queue:
        board, path = queue.popleft()
        
        # Check if solved
        if not any(c == '@' for c in board):
            return path
        
        player_pos = get_pos(board, '*%')[0]
        
        for move, (dx, dy) in moves.items():
            valid, box_pos = is_valid_move(board, player_pos, dx, dy)
            if valid:
                new_board = make_move(board, player_pos, dx, dy)
                if new_board not in visited:
                    visited.add(new_board)
                    queue.append((new_board, path + move))
    
    return None

# Initialize the board
board = """+++++++++
+-X---@%+
+-----@@+
+--@--+X+
+++X$-+++
+++++++++"""

solution = solve_sokoban(board)
print(solution)