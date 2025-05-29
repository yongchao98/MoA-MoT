from collections import deque

# Define the initial state of the puzzle
initial_state = [
    "+ + + + + +",
    "+ + - * - +",
    "+ $ $ @ @ +",
    "+ X - - X +",
    "+ X @ - - +",
    "+ - - - - +",
    "+ - - - - +",
    "+ - - - - +",
    "+ - - @ X +",
    "+ + + + + +"
]

# Convert the state to a more manageable format
state = [list(row.replace(" ", "")) for row in initial_state]

# Define the directions for movement
directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}

# Find the initial position of the player
def find_player(state):
    for i, row in enumerate(state):
        for j, cell in enumerate(row):
            if cell == '*' or cell == '%':
                return i, j
    return None

# Check if the move is valid
def is_valid_move(state, x, y, dx, dy):
    if state[x + dx][y + dy] in '-X':  # Move to empty or goal
        return True
    if state[x + dx][y + dy] in '@$':  # Move to box
        if state[x + 2*dx][y + 2*dy] in '-X':  # Check if box can be pushed
            return True
    return False

# Apply the move to the state
def apply_move(state, x, y, dx, dy):
    new_state = [row[:] for row in state]
    if new_state[x + dx][y + dy] in '-X':  # Move to empty or goal
        new_state[x][y], new_state[x + dx][y + dy] = '-', '*'
    elif new_state[x + dx][y + dy] in '@$':  # Move to box
        new_state[x][y], new_state[x + dx][y + dy], new_state[x + 2*dx][y + 2*dy] = '-', '*', '@'
    return new_state

# Check if the puzzle is solved
def is_solved(state):
    for row in state:
        if '@' in row:
            return False
    return True

# Perform a breadth-first search to find the solution
def solve_sokoban(state):
    start = find_player(state)
    queue = deque([(state, start, "")])
    visited = set()
    visited.add(tuple(tuple(row) for row in state))
    
    while queue:
        current_state, (x, y), path = queue.popleft()
        
        if is_solved(current_state):
            return path
        
        for move, (dx, dy) in directions.items():
            if is_valid_move(current_state, x, y, dx, dy):
                new_state = apply_move(current_state, x, y, dx, dy)
                new_tuple_state = tuple(tuple(row) for row in new_state)
                if new_tuple_state not in visited:
                    visited.add(new_tuple_state)
                    queue.append((new_state, (x + dx, y + dy), path + move))
    
    return "No solution"

# Find the solution
solution = solve_sokoban(state)
print(solution)