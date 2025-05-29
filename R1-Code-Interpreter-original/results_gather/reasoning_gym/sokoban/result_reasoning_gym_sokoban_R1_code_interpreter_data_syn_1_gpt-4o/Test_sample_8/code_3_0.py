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

# Find the initial position of the player and boxes
def find_positions(board):
    player = None
    boxes = []
    goals = []
    for i, row in enumerate(board):
        for j, cell in enumerate(row):
            if cell == '*' or cell == '%':
                player = (i, j)
            if cell == '@' or cell == '$':
                boxes.append((i, j))
            if cell == 'X' or cell == '$':
                goals.append((i, j))
    return player, boxes, goals

# Check if a move is valid
def is_valid_move(board, player, boxes, dx, dy):
    x, y = player
    new_x, new_y = x + dx, y + dy
    if board[new_x][new_y] in '-X':
        return True
    if (new_x, new_y) in boxes:
        new_box_x, new_box_y = new_x + dx, new_y + dy
        if board[new_box_x][new_box_y] in '-X' and (new_box_x, new_box_y) not in boxes:
            return True
    return False

# Apply a move to the board
def apply_move(player, boxes, dx, dy):
    x, y = player
    new_x, new_y = x + dx, y + dy
    new_boxes = set(boxes)
    if (new_x, new_y) in new_boxes:
        new_boxes.remove((new_x, new_y))
        new_boxes.add((new_x + dx, new_y + dy))
    return (new_x, new_y), frozenset(new_boxes)

# Check if the board is in a solved state
def is_solved(boxes, goals):
    return all(box in goals for box in boxes)

# Heuristic function: sum of Manhattan distances from each box to the nearest goal
def heuristic(boxes, goals):
    total_distance = 0
    for box in boxes:
        min_distance = float('inf')
        for goal in goals:
            distance = abs(box[0] - goal[0]) + abs(box[1] - goal[1])
            min_distance = min(min_distance, distance)
        total_distance += min_distance
    return total_distance

# Perform A* search to find the solution
def solve_sokoban(board):
    player, boxes, goals = find_positions(board)
    queue = []
    heappush(queue, (0, player, frozenset(boxes), ""))
    visited = set()
    visited.add((player, frozenset(boxes)))

    while queue:
        _, current_player, current_boxes, path = heappop(queue)

        if is_solved(current_boxes, goals):
            return path

        for move, (dx, dy) in directions.items():
            if is_valid_move(board, current_player, current_boxes, dx, dy):
                new_player, new_boxes = apply_move(current_player, current_boxes, dx, dy)
                new_state = (new_player, new_boxes)
                if new_state not in visited:
                    visited.add(new_state)
                    cost = len(path) + 1 + heuristic(new_boxes, goals)
                    heappush(queue, (cost, new_player, new_boxes, path + move))

    return "No solution"

# Find the solution
solution = solve_sokoban(board)
print(solution)