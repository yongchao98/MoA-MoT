from heapq import heappop, heappush

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

# Find all goal positions
def find_goals(board):
    goals = []
    for i, row in enumerate(board):
        for j, cell in enumerate(row):
            if cell == 'X' or cell == '$':
                goals.append((i, j))
    return goals

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

# Heuristic function: sum of Manhattan distances from each box to the nearest goal
def heuristic(board, goals):
    total_distance = 0
    for i, row in enumerate(board):
        for j, cell in enumerate(row):
            if cell == '@':
                min_distance = float('inf')
                for gx, gy in goals:
                    distance = abs(i - gx) + abs(j - gy)
                    min_distance = min(min_distance, distance)
                total_distance += min_distance
    return total_distance

# Check for deadlocks
def is_deadlock(board):
    # Simple deadlock detection: check for boxes in corners not on goals
    for i, row in enumerate(board):
        for j, cell in enumerate(row):
            if cell == '@':
                if (board[i-1][j] == '+' or board[i+1][j] == '+') and (board[i][j-1] == '+' or board[i][j+1] == '+'):
                    return True
    return False

# Perform A* search to find the solution
def solve_sokoban(board):
    start = find_player(board)
    goals = find_goals(board)
    queue = []
    heappush(queue, (0, board, start, ""))
    visited = set()
    visited.add(tuple(map(tuple, board)))

    while queue:
        _, current_board, (x, y), path = heappop(queue)

        if is_solved(current_board):
            return path

        for move, (dx, dy) in directions.items():
            if is_valid_move(current_board, x, y, dx, dy):
                new_board = apply_move(current_board, x, y, dx, dy)
                new_state = tuple(map(tuple, new_board))
                if new_state not in visited and not is_deadlock(new_board):
                    visited.add(new_state)
                    cost = len(path) + 1 + heuristic(new_board, goals)
                    heappush(queue, (cost, new_board, (x + dx, y + dy), path + move))

    return "No solution"

# Find the solution
solution = solve_sokoban(board)
print(solution)