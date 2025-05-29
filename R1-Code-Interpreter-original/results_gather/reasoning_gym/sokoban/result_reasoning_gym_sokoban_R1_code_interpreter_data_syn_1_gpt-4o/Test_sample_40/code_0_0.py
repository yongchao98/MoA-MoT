from collections import deque

# Define the initial state
initial_state = {
    'player': (1, 1),
    'boxes': {(1, 2), (2, 1), (2, 2)},
    'goals': {(1, 3), (3, 1), (6, 3)}
}

# Define the grid
grid = [
    "+ + + + + +",
    "+ * @ - X +",
    "+ @ @ $ + +",
    "+ X - $ + +",
    "+ + - + $ +",
    "+ - - + + +",
    "+ - - X + +",
    "+ - + + + +",
    "+ + + + + +",
    "+ + + + + +"
]

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

# Check if a position is within the grid and not a wall
def is_valid_position(pos):
    x, y = pos
    return 0 <= x < len(grid) and 0 <= y < len(grid[0]) and grid[x][y] != '+'

# Check if all boxes are on goals
def is_solved(boxes, goals):
    return boxes == goals

# Perform BFS to find the solution
def solve_sokoban(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()

    while queue:
        state, path = queue.popleft()
        player = state['player']
        boxes = state['boxes']

        if is_solved(boxes, initial_state['goals']):
            return path

        for move, (dx, dy) in moves.items():
            new_player = (player[0] + dx, player[1] + dy)

            if not is_valid_position(new_player):
                continue

            if new_player in boxes:
                new_box = (new_player[0] + dx, new_player[1] + dy)
                if not is_valid_position(new_box) or new_box in boxes:
                    continue
                new_boxes = boxes - {new_player} | {new_box}
            else:
                new_boxes = boxes

            new_state = {'player': new_player, 'boxes': new_boxes}

            if (new_player, frozenset(new_boxes)) not in visited:
                visited.add((new_player, frozenset(new_boxes)))
                queue.append((new_state, path + move))

    return None

solution = solve_sokoban(initial_state)
print(solution)