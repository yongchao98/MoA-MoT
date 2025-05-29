from collections import deque

# Define the initial state of the board
initial_board = [
    "++++++++",
    "+*@X---+",
    "+@@$---+",
    "+X-$+++",
    "++X++++",
    "++++++++",
    "++++++++",
    "+++$+++",
    "++++++++"
]

# Define the goal state
goal_positions = {(1, 3), (3, 1), (3, 2)}

# Define the possible moves
moves = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}

# Find the initial position of the player
def find_player(board):
    for r, row in enumerate(board):
        for c, cell in enumerate(row):
            if cell == '*' or cell == '%':
                return r, c
    return None

# Check if the current state is a goal state
def is_goal_state(board):
    for r, row in enumerate(board):
        for c, cell in enumerate(row):
            if cell == '@' and (r, c) not in goal_positions:
                return False
    return True

# Move the player and update the board
def move_player(board, player_pos, direction):
    r, c = player_pos
    dr, dc = moves[direction]
    new_r, new_c = r + dr, c + dc
    if board[new_r][new_c] in '-X':  # Move to empty or goal
        new_board = [list(row) for row in board]
        new_board[r][c] = '-' if board[r][c] == '*' else 'X'
        new_board[new_r][new_c] = '*' if board[new_r][new_c] == '-' else '%'
        return new_board, (new_r, new_c)
    elif board[new_r][new_c] in '@$':  # Push a box
        box_r, box_c = new_r + dr, new_c + dc
        if board[box_r][box_c] in '-X':  # Box can be pushed
            new_board = [list(row) for row in board]
            new_board[r][c] = '-' if board[r][c] == '*' else 'X'
            new_board[new_r][new_c] = '*' if board[new_r][new_c] == '@' else '%'
            new_board[box_r][box_c] = '@' if board[box_r][box_c] == '-' else '$'
            return new_board, (new_r, new_c)
    return None, player_pos

# Perform BFS to find the solution
def solve_sokoban(initial_board):
    initial_player_pos = find_player(initial_board)
    queue = deque([(initial_board, initial_player_pos, "")])
    visited = set()
    visited.add((tuple(initial_board), initial_player_pos))
    
    while queue:
        board, player_pos, path = queue.popleft()
        
        if is_goal_state(board):
            return path
        
        for direction in moves:
            new_board, new_player_pos = move_player(board, player_pos, direction)
            if new_board and (tuple(new_board), new_player_pos) not in visited:
                visited.add((tuple(new_board), new_player_pos))
                queue.append((new_board, new_player_pos, path + direction))
    
    return None

# Solve the puzzle
solution = solve_sokoban(initial_board)
print(solution)