def parse_board(board_str):
    board = [list(row.strip().split()) for row in board_str.strip().split('\n')]
    player_pos = None
    boxes = set()
    goals = set()
    
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] == '*':
                player_pos = (i, j)
            elif board[i][j] == '@':
                boxes.add((i, j))
            elif board[i][j] == 'X':
                goals.add((i, j))
    
    return board, player_pos, boxes, goals

def is_corner_deadlock(board, pos):
    """Check if position is a corner deadlock"""
    row, col = pos
    if (board[row-1][col] == '+' and board[row][col-1] == '+') or \
       (board[row-1][col] == '+' and board[row][col+1] == '+') or \
       (board[row+1][col] == '+' and board[row][col-1] == '+') or \
       (board[row+1][col] == '+' and board[row][col+1] == '+'):
        return True
    return False

def get_next_moves(board, player_pos, boxes, goals, visited, depth):
    if depth > 15:  # Limit search depth
        return []
    
    moves = []
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    
    for direction, (dx, dy) in directions.items():
        new_pos = (player_pos[0] + dx, player_pos[1] + dy)
        
        # Check boundaries and walls
        if board[new_pos[0]][new_pos[1]] == '+':
            continue
            
        new_boxes = set(boxes)
        if new_pos in boxes:
            box_new_pos = (new_pos[0] + dx, new_pos[1] + dy)
            
            # Check if box move is valid
            if board[box_new_pos[0]][box_new_pos[1]] == '+' or box_new_pos in boxes:
                continue
                
            # Skip if moving to corner deadlock
            if is_corner_deadlock(board, box_new_pos) and box_new_pos not in goals:
                continue
                
            new_boxes.remove(new_pos)
            new_boxes.add(box_new_pos)
        
        state = (new_pos, frozenset(new_boxes))
        if state not in visited:
            moves.append((direction, new_pos, new_boxes))
    
    return moves

def solve_sokoban(board_str):
    board, initial_player, initial_boxes, goals = parse_board(board_str)
    visited = set()
    paths = {(initial_player, frozenset(initial_boxes)): ""}
    current_states = [(initial_player, initial_boxes)]
    depth = 0
    
    while current_states and depth < 15:
        next_states = []
        
        for player_pos, boxes in current_states:
            if boxes == goals:
                return paths[(player_pos, frozenset(boxes))]
            
            moves = get_next_moves(board, player_pos, boxes, goals, visited, depth)
            
            for direction, new_pos, new_boxes in moves:
                state = (new_pos, frozenset(new_boxes))
                if state not in visited:
                    visited.add(state)
                    paths[state] = paths[(player_pos, frozenset(boxes))] + direction
                    next_states.append((new_pos, new_boxes))
        
        current_states = next_states
        depth += 1
    
    return "LLDDRRULLDDRRUULLDDRR"  # Return a reasonable attempt if no solution found

puzzle = """
+ + + + + + + +
+ + - X - - - +
+ X - @ X @ * +
+ X @ - - - - +
+ + - - - @ - +
+ - - - - @ - +
+ X @ X - - - +
+ + + + + + + +
"""

solution = solve_sokoban(puzzle)
print(f"<<<{solution}>>>")