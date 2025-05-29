from collections import deque

# Define the initial state of the puzzle
initial_state = [
    ['+', '+', '+', '+', '+', '+'],
    ['+', 'X', '-', '-', '-', '+'],
    ['+', '-', '@', '-', '$', '+'],
    ['+', '-', '$', '+', '+', '+'],
    ['+', '-', '@', 'X', '-', '+'],
    ['+', '-', '*', '$', '-', '+'],
    ['+', '-', '-', '-', '-', '+'],
    ['+', '+', '+', '+', '+', '+']
]

# Define the possible moves
moves = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}

# Find the player's initial position
def find_player(state):
    for i, row in enumerate(state):
        for j, cell in enumerate(row):
            if cell == '*':
                return i, j
    return None

# Check if the current state is the goal state
def is_goal(state):
    for i in range(len(state)):
        for j in range(len(state[i])):
            if state[i][j] == '@':
                return False
    return True

# Apply a move to the current state
def apply_move(state, move):
    new_state = [row[:] for row in state]
    pi, pj = find_player(state)
    di, dj = moves[move]
    ni, nj = pi + di, pj + dj
    if new_state[ni][nj] in ('-', 'X'):  # Move to empty or goal
        new_state[pi][pj], new_state[ni][nj] = '-', '*'
    elif new_state[ni][nj] in ('@', '$'):  # Push a box
        bi, bj = ni + di, nj + dj
        if new_state[bi][bj] in ('-', 'X'):  # Box can be pushed
            new_state[pi][pj], new_state[ni][nj], new_state[bi][bj] = '-', '*', '@' if new_state[bi][bj] == '-' else '$'
    return new_state

# Perform a breadth-first search to find the solution
def solve_sokoban():
    queue = deque([(initial_state, "")])
    visited = set()
    while queue:
        state, path = queue.popleft()
        if is_goal(state):
            return path
        for move in moves:
            new_state = apply_move(state, move)
            state_tuple = tuple(tuple(row) for row in new_state)
            if state_tuple not in visited:
                visited.add(state_tuple)
                queue.append((new_state, path + move))
    return None

# Find the solution
solution = solve_sokoban()
print(solution)