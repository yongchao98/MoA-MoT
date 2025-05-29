from collections import deque

# Define the initial state of the board
initial_board = [
    "+ + + + + + + + +",
    "+ * $ - - - + + +",
    "+ @ - - - - + + +",
    "+ X @ @ - - - + +",
    "+ + X $ - X - + +",
    "+ + + X - @ @ + +",
    "+ $ + X - - - - +",
    "+ + - @ - - - - +",
    "+ + - - - - - X +",
    "+ + + + + + + + +"
]

# Convert the board to a more manageable format
board = [list(row.replace(" ", "")) for row in initial_board]

# Define the directions for movement
directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}

# Find the initial position of the player
def find_player(board):
    for i, row in enumerate(board):
        for j, cell in enumerate(row):
            if cell == '*' or cell == '%':
                return i, j
    return None

# Check if a move is valid
def is_valid_move(board, x, y, dx, dy):
    if board[x + dx][y + dy] in '-X':
        return True
    if board[x + dx][y + dy] in '@$':
        if board[x + 2*dx][y + 2*dy] in '-X':
            return True
    return False

# Apply a move to the board
def apply_move(board, x, y, dx, dy):
    new_board = [row[:] for row in board]
    if new_board[x + dx][y + dy] in '-X':
        new_board[x][y], new_board[x + dx][y + dy] = '-', '*'
    elif new_board[x + dx][y + dy] in '@$':
        new_board[x][y], new_board[x + dx][y + dy], new_board[x + 2*dx][y + 2*dy] = '-', '*', '@'
    return new_board

# Check if the board is in a solved state
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
        current_board, (x, y), path = queue.popleft()

        if is_solved(current_board):
            return path

        for move, (dx, dy) in directions.items():
            if is_valid_move(current_board, x, y, dx, dy):
                new_board = apply_move(current_board, x, y, dx, dy)
                new_state = tuple(map(tuple, new_board))
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_board, (x + dx, y + dy), path + move))

    return "No solution"

# Find the solution
solution = solve_sokoban(board)
print(solution)