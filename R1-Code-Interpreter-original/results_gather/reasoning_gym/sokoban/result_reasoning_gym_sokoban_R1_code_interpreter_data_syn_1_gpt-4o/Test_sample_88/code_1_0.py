from heapq import heappop, heappush

# Define the initial state of the puzzle
initial_state = [
    "+ + + + + + + + + +",
    "+ - - - - - - - - +",
    "+ - X - $ $ - @ - +",
    "+ - - - X X X - - +",
    "+ + X - @ @ @ @ * +",
    "+ + - - - - - - - +",
    "+ + + + + + + + + +"
]

# Convert the state to a more manageable format
def parse_state(state):
    return [list(row.replace(" ", "")) for row in state]

# Find the player's position
def find_player(state):
    for y, row in enumerate(state):
        for x, cell in enumerate(row):
            if cell == '*':
                return (x, y)
    return None

# Check if a move is valid
def is_valid_move(state, x, y, dx, dy):
    target_x, target_y = x + dx, y + dy
    if state[target_y][target_x] in ('-', 'X'):
        return True
    if state[target_y][target_x] in ('@', '$'):
        beyond_x, beyond_y = target_x + dx, target_y + dy
        if state[beyond_y][beyond_x] in ('-', 'X'):
            return True
    return False

# Apply a move to the state
def apply_move(state, x, y, dx, dy):
    new_state = [row[:] for row in state]
    target_x, target_y = x + dx, y + dy
    if new_state[target_y][target_x] in ('-', 'X'):
        new_state[y][x], new_state[target_y][target_x] = '-', '*'
    elif new_state[target_y][target_x] in ('@', '$'):
        beyond_x, beyond_y = target_x + dx, target_y + dy
        if new_state[beyond_y][beyond_x] in ('-', 'X'):
            new_state[y][x], new_state[target_y][target_x] = '-', '*'
            new_state[beyond_y][beyond_x] = '@' if new_state[beyond_y][beyond_x] == '-' else '$'
    return new_state

# Check if the puzzle is solved
def is_solved(state):
    for row in state:
        for cell in row:
            if cell == '@':
                return False
    return True

# Heuristic function for A* search
def heuristic(state):
    # Count the number of boxes not on goals
    return sum(row.count('@') for row in state)

# Perform an A* search to find the solution
def solve_sokoban(initial_state):
    state = parse_state(initial_state)
    player_x, player_y = find_player(state)
    queue = [(0, state, player_x, player_y, "")]
    visited = set()

    while queue:
        _, current_state, x, y, path = heappop(queue)
        if is_solved(current_state):
            return path

        for dx, dy, move in [(-1, 0, 'L'), (1, 0, 'R'), (0, -1, 'U'), (0, 1, 'D')]:
            if is_valid_move(current_state, x, y, dx, dy):
                new_state = apply_move(current_state, x, y, dx, dy)
                new_player_x, new_player_y = x + dx, y + dy
                state_tuple = tuple(tuple(row) for row in new_state)
                if state_tuple not in visited:
                    visited.add(state_tuple)
                    cost = len(path) + 1 + heuristic(new_state)
                    heappush(queue, (cost, new_state, new_player_x, new_player_y, path + move))

    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)