from heapq import heappop, heappush

# Define the initial state of the puzzle
initial_state = [
    "++++++++++",
    "++-----+$+",
    "+-------++",
    "+X-X---$++",
    "+-@----+++",
    "+X-@@@++++",
    "++X-@%$+++",
    "++++++++++"
]

# Define the goal state
goal_positions = {(3, 1), (3, 3), (5, 1), (6, 2)}

# Define the possible moves
moves = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}

# Find the player's initial position
def find_player(state):
    for i, row in enumerate(state):
        for j, cell in enumerate(row):
            if cell == '%' or cell == '*':
                return i, j
    return None

# Check if the current state is the goal state
def is_goal(state):
    for i, row in enumerate(state):
        for j, cell in enumerate(row):
            if cell == '@' and (i, j) not in goal_positions:
                return False
    return True

# Generate the next state after a move
def move(state, direction):
    new_state = [list(row) for row in state]
    pi, pj = find_player(state)
    di, dj = moves[direction]
    ni, nj = pi + di, pj + dj
    if new_state[ni][nj] in '-X':  # Move to empty or goal
        new_state[pi][pj] = '-' if new_state[pi][pj] == '*' else 'X'
        new_state[ni][nj] = '%' if new_state[ni][nj] == 'X' else '*'
    elif new_state[ni][nj] in '@$':  # Push box
        nni, nnj = ni + di, nj + dj
        if new_state[nni][nnj] in '-X':  # Box can be pushed
            new_state[pi][pj] = '-' if new_state[pi][pj] == '*' else 'X'
            new_state[ni][nj] = '%' if new_state[ni][nj] == 'X' else '*'
            new_state[nni][nnj] = '$' if new_state[nni][nnj] == 'X' else '@'
    return [''.join(row) for row in new_state]

# Heuristic function: sum of Manhattan distances from each box to the nearest goal
def heuristic(state):
    boxes = [(i, j) for i, row in enumerate(state) for j, cell in enumerate(row) if cell in '@$']
    return sum(min(abs(bi - gi) + abs(bj - gj) for gi, gj in goal_positions) for bi, bj in boxes)

# Perform an A* search to find the solution
def solve_sokoban():
    queue = []
    heappush(queue, (0, initial_state, ""))
    visited = set()
    while queue:
        _, state, path = heappop(queue)
        if is_goal(state):
            return path
        if tuple(state) in visited:
            continue
        visited.add(tuple(state))
        for direction in moves:
            new_state = move(state, direction)
            if tuple(new_state) not in visited:
                cost = len(path) + 1 + heuristic(new_state)
                heappush(queue, (cost, new_state, path + direction))
    return None

# Find the solution
solution = solve_sokoban()
print(solution)