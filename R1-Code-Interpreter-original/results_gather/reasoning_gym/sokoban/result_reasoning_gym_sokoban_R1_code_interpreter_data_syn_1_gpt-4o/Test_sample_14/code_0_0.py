from collections import deque

# Define the initial state of the puzzle
initial_state = [
    ['+', '+', '+', '+', '+', '+'],
    ['+', 'X', '-', '-', '-', '+'],
    ['+', '+', '-', '-', '-', '+'],
    ['+', '-', '-', '@', '-', '+'],
    ['+', '-', '-', '-', '-', '+'],
    ['+', '-', 'X', '@', '-', '+'],
    ['+', '$', '$', '-', '-', '+'],
    ['+', '*', '@', 'X', '-', '+'],
    ['+', '+', '+', '+', '+', '+']
]

# Define the goal state
goal_state = [
    ['+', '+', '+', '+', '+', '+'],
    ['+', 'X', '-', '-', '-', '+'],
    ['+', '+', '-', '-', '-', '+'],
    ['+', '-', '-', '-', '-', '+'],
    ['+', '-', '-', '-', '-', '+'],
    ['+', '-', 'X', '-', '-', '+'],
    ['+', '$', '$', '-', '-', '+'],
    ['+', '-', '$', 'X', '-', '+'],
    ['+', '+', '+', '+', '+', '+']
]

# Define the possible moves
moves = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}

# Find the player's initial position
def find_player(state):
    for r, row in enumerate(state):
        for c, cell in enumerate(row):
            if cell == '*' or cell == '%':
                return r, c
    return None

# Check if the current state is the goal state
def is_goal(state):
    for r in range(len(state)):
        for c in range(len(state[0])):
            if state[r][c] == '@' and goal_state[r][c] != 'X':
                return False
    return True

# Apply a move to the current state
def apply_move(state, move):
    r, c = find_player(state)
    dr, dc = moves[move]
    nr, nc = r + dr, c + dc
    nnr, nnc = nr + dr, nc + dc

    # Check if the move is valid
    if state[nr][nc] in ('-', 'X'):
        # Move the player
        new_state = [row[:] for row in state]
        new_state[r][c] = '-' if state[r][c] == '*' else 'X'
        new_state[nr][nc] = '*' if state[nr][nc] == '-' else '%'
        return new_state
    elif state[nr][nc] in ('@', '$') and state[nnr][nnc] in ('-', 'X'):
        # Push the box
        new_state = [row[:] for row in state]
        new_state[r][c] = '-' if state[r][c] == '*' else 'X'
        new_state[nr][nc] = '*' if state[nr][nc] == '@' else '%'
        new_state[nnr][nnc] = '@' if state[nnr][nnc] == '-' else '$'
        return new_state
    return None

# Perform a breadth-first search to find the solution
def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()
    visited.add(tuple(map(tuple, initial_state)))

    while queue:
        current_state, path = queue.popleft()

        if is_goal(current_state):
            return path

        for move in moves:
            new_state = apply_move(current_state, move)
            if new_state and tuple(map(tuple, new_state)) not in visited:
                visited.add(tuple(map(tuple, new_state)))
                queue.append((new_state, path + move))

    return None

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)