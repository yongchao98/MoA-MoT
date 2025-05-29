from heapq import heappop, heappush

# Define the initial state of the board
initial_board = [
    "+ + + + + + + + + +",
    "+ - * - - - - - X +",
    "+ @ @ @ @ - - @ - +",
    "+ - - - - X - $ - +",
    "+ X - - - $ X - - +",
    "+ X - - - - - - - +",
    "+ + + + + + + + + +"
]

# Convert the board into a more manageable format
board = [list(row.split()) for row in initial_board]

# Define the directions for movement
directions = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

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

# Find all box positions
def find_boxes(board):
    boxes = []
    for i, row in enumerate(board):
        for j, cell in enumerate(row):
            if cell == '@' or cell == '$':
                boxes.append((i, j))
    return boxes

# Check if a move is valid
def is_valid_move(board, x, y, dx, dy):
    if board[x + dx][y + dy] in ('-', 'X'):  # Move to empty or goal
        return True
    if board[x + dx][y + dy] in ('@', '$'):  # Move to box
        if board[x + 2*dx][y + 2*dy] in ('-', 'X'):  # Check if box can be pushed
            return True
    return False

# Apply a move to the board
def apply_move(board, x, y, dx, dy):
    new_board = [row[:] for row in board]
    if new_board[x + dx][y + dy] in ('-', 'X'):  # Move to empty or goal
        new_board[x][y] = '-' if new_board[x][y] == '*' else 'X'
        new_board[x + dx][y + dy] = '*' if new_board[x + dx][y + dy] == '-' else '%'
    elif new_board[x + dx][y + dy] in ('@', '$'):  # Move to box
        new_board[x][y] = '-' if new_board[x][y] == '*' else 'X'
        new_board[x + dx][y + dy] = '*' if new_board[x + dx][y + dy] == '@' else '%'
        new_board[x + 2*dx][y + 2*dy] = '@' if new_board[x + 2*dx][y + 2*dy] == '-' else '$'
    return new_board

# Check if the board is in a winning state
def is_winning(board):
    for row in board:
        for cell in row:
            if cell == '@':  # Any box not on a goal
                return False
    return True

# Calculate the heuristic: sum of Manhattan distances from boxes to nearest goals
def heuristic(boxes, goals):
    total_distance = 0
    for bx, by in boxes:
        min_distance = float('inf')
        for gx, gy in goals:
            distance = abs(bx - gx) + abs(by - gy)
            min_distance = min(min_distance, distance)
        total_distance += min_distance
    return total_distance

# Perform A* search to find the solution
def solve_sokoban(board):
    start = find_player(board)
    goals = find_goals(board)
    boxes = find_boxes(board)
    queue = []
    heappush(queue, (0, board, start, "", boxes))
    visited = set()
    visited.add((tuple(map(tuple, board)), tuple(boxes)))

    while queue:
        _, current_board, (x, y), path, current_boxes = heappop(queue)

        if is_winning(current_board):
            return path

        for move, (dx, dy) in directions.items():
            if is_valid_move(current_board, x, y, dx, dy):
                new_board = apply_move(current_board, x, y, dx, dy)
                new_boxes = find_boxes(new_board)
                new_state = (tuple(map(tuple, new_board)), tuple(new_boxes))
                if new_state not in visited:
                    visited.add(new_state)
                    new_cost = len(path) + 1 + heuristic(new_boxes, goals)
                    heappush(queue, (new_cost, new_board, (x + dx, y + dy), path + move, new_boxes))

    return "No solution"

# Find the solution
solution = solve_sokoban(board)
print(solution)