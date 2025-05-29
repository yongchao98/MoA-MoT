from collections import deque

def is_valid(board, row, col):
    return 0 <= row < len(board) and 0 <= col < len(board[0]) and board[row][col] != '+'

def get_next_states(board, player_pos):
    moves = [('U', -1, 0), ('D', 1, 0), ('L', 0, -1), ('R', 0, 1)]
    next_states = []
    
    for move, dr, dc in moves:
        r, c = player_pos
        new_r, new_c = r + dr, c + dc
        
        if not is_valid(board, new_r, new_c):
            continue
            
        new_board = [row[:] for row in board]
        
        # Moving to empty space or goal
        if new_board[new_r][new_c] in '-X':
            new_board[r][c] = '-' if new_board[r][c] == '*' else 'X'
            new_board[new_r][new_c] = '*' if new_board[new_r][new_c] == '-' else '%'
            next_states.append((move, new_board, (new_r, new_c)))
            
        # Moving a box
        elif new_board[new_r][new_c] in '@$':
            box_r, box_c = new_r + dr, new_c + dc
            if is_valid(board, box_r, box_c) and new_board[box_r][box_c] in '-X':
                new_board[r][c] = '-' if new_board[r][c] == '*' else 'X'
                new_board[new_r][new_c] = '*' if new_board[new_r][new_c] == '@' else '%'
                new_board[box_r][box_c] = '$' if new_board[box_r][box_c] == 'X' else '@'
                next_states.append((move, new_board, (new_r, new_c)))
    
    return next_states

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def count_goals_and_boxes_on_goals(board):
    goals = 0
    boxes_on_goals = 0
    for row in board:
        for cell in row:
            if cell in 'X%':
                goals += 1
            if cell in '$%':
                boxes_on_goals += 1
    return goals, boxes_on_goals

def is_goal(board):
    goals, boxes_on_goals = count_goals_and_boxes_on_goals(board)
    return goals == 0 or (goals > 0 and goals == boxes_on_goals)

def solve_sokoban(initial_board):
    # Parse the board
    initial_state = [list(line.strip().replace(' ', '')) for line in initial_board.strip().split('\n')]
    
    # Find player position
    player_pos = None
    for i in range(len(initial_state)):
        for j in range(len(initial_state[i])):
            if initial_state[i][j] in '*%':
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    visited = set()
    queue = deque([([], initial_state, player_pos)])
    max_iterations = 10000  # Prevent infinite loops
    iterations = 0
    
    while queue and iterations < max_iterations:
        path, current_board, pos = queue.popleft()
        board_str = board_to_string(current_board)
        iterations += 1
        
        if board_str in visited:
            continue
            
        visited.add(board_str)
        
        if is_goal(current_board):
            return ''.join(path)
            
        for move, new_board, new_pos in get_next_states(current_board, pos):
            queue.append((path + [move], new_board, new_pos))
    
    return None

# Initial board
initial_board = """++++++++
+XX@-@-+
+--XX@-+
+--$@-+
+-*---+
++++++++"""

solution = solve_sokoban(initial_board)
print(solution)