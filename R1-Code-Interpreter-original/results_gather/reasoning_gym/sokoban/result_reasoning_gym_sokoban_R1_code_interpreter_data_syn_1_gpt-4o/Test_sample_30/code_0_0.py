from collections import deque

# Define the initial state of the board
initial_board = [
    ['+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '%', '@', 'X', '-', '-', '-', '+'],
    ['+', '@', '-', '@', '@', '-', '-', '+'],
    ['+', '-', '-', '-', '-', '-', '-', '+'],
    ['+', '-', '-', 'X', 'X', '-', '-', '+'],
    ['+', '-', '-', '-', '@', 'X', '-', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+']
]

# Define the goal state
goal_positions = {(1, 3), (4, 3), (4, 4), (5, 5)}

# Define the directions for movement
directions = {'L': (0, -1), 'R': (0, 1), 'U': (-1, 0), 'D': (1, 0)}

# Find the initial position of the player
def find_player(board):
    for r, row in enumerate(board):
        for c, cell in enumerate(row):
            if cell == '%' or cell == '-':
                return (r, c)
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
    dr, dc = directions[direction]
    new_r, new_c = r + dr, c + dc
    if board[new_r][new_c] in ('-', 'X'):  # Move to empty or goal
        board[r][c], board[new_r][new_c] = '-', '%'
        return (new_r, new_c), board
    elif board[new_r][new_c] == '@':  # Push box
        box_new_r, box_new_c = new_r + dr, new_c + dc
        if board[box_new_r][box_new_c] in ('-', 'X'):
            board[r][c], board[new_r][new_c], board[box_new_r][box_new_c] = '-', '%', '@'
            return (new_r, new_c), board
    return None, board

# Perform a breadth-first search to find the solution
def solve_sokoban(initial_board):
    initial_player_pos = find_player(initial_board)
    queue = deque([(initial_player_pos, initial_board, "")])
    visited = set()
    visited.add((initial_player_pos, tuple(map(tuple, initial_board))))

    while queue:
        player_pos, board, path = queue.popleft()
        if is_goal_state(board):
            return path
        for direction in directions:
            new_player_pos, new_board = move_player([row[:] for row in board], player_pos, direction)
            if new_player_pos and (new_player_pos, tuple(map(tuple, new_board))) not in visited:
                visited.add((new_player_pos, tuple(map(tuple, new_board))))
                queue.append((new_player_pos, new_board, path + direction))
    return None

# Find the solution
solution = solve_sokoban(initial_board)
print(solution)