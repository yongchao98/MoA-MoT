from collections import deque

# Define the initial state of the board
initial_board = [
    "+ + + + + + +",
    "+ $ + - - X +",
    "+ - - - - - +",
    "+ - - - @ - +",
    "+ - X - @ - +",
    "+ - - - X - +",
    "+ - - @ X X +",
    "+ - - - @ @ +",
    "+ - - - @ % +",
    "+ + + + + + +"
]

# Convert the board into a more manageable format
board = [list(row.split()) for row in initial_board]

# Define the directions for movement
directions = {'L': (0, -1), 'R': (0, 1), 'U': (-1, 0), 'D': (1, 0)}

# Find the initial position of the player
def find_player(board):
    for r, row in enumerate(board):
        for c, cell in enumerate(row):
            if cell == '*' or cell == '%':
                return r, c
    return None

# Check if a move is valid
def is_valid_move(board, r, c, dr, dc):
    if board[r + dr][c + dc] in ('-', 'X'):  # Move to empty or goal
        return True
    if board[r + dr][c + dc] in ('@', '$'):  # Move to box
        if board[r + 2*dr][c + 2*dc] in ('-', 'X'):  # Box can be pushed
            return True
    return False

# Apply a move to the board
def apply_move(board, r, c, dr, dc):
    new_board = [row[:] for row in board]
    if new_board[r + dr][c + dc] in ('-', 'X'):  # Move to empty or goal
        new_board[r][c] = '-' if new_board[r][c] == '*' else 'X'
        new_board[r + dr][c + dc] = '*' if new_board[r + dr][c + dc] == '-' else '%'
    elif new_board[r + dr][c + dc] in ('@', '$'):  # Move to box
        new_board[r][c] = '-' if new_board[r][c] == '*' else 'X'
        new_board[r + dr][c + dc] = '*' if new_board[r + dr][c + dc] == '@' else '%'
        new_board[r + 2*dr][c + 2*dc] = '@' if new_board[r + 2*dr][c + 2*dc] == '-' else '$'
    return new_board

# Check if the board is solved
def is_solved(board):
    for row in board:
        for cell in row:
            if cell == '@':
                return False
    return True

# Perform BFS to find the solution
def solve_sokoban(board):
    start = find_player(board)
    queue = deque([(board, start, "")])
    visited = set()
    visited.add(tuple(map(tuple, board)))
    
    while queue:
        current_board, (r, c), path = queue.popleft()
        
        if is_solved(current_board):
            return path
        
        for move, (dr, dc) in directions.items():
            if is_valid_move(current_board, r, c, dr, dc):
                new_board = apply_move(current_board, r, c, dr, dc)
                new_state = tuple(map(tuple, new_board))
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_board, (r + dr, c + dc), path + move))
    
    return "No solution"

# Find the solution
solution = solve_sokoban(board)
print(solution)