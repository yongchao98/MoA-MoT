from collections import deque

# Define the initial state of the puzzle
initial_state = [
    "+ + + + + + + + +",
    "+ * @ - - - X $ +",
    "+ @ - @ - X + + +",
    "+ - - - - $ $ $ +",
    "+ - + + + + + + +",
    "+ - - + - - - - +",
    "+ X @ - - + + + +",
    "+ + - @ X + + + +",
    "+ + X + + + + + +",
    "+ + + + + + + + +"
]

# Convert the initial state into a more manageable format
def parse_state(state):
    grid = [list(row.split()) for row in state]
    player_pos = None
    boxes = set()
    goals = set()
    
    for r, row in enumerate(grid):
        for c, cell in enumerate(row):
            if cell == '*':
                player_pos = (r, c)
            elif cell == '@':
                boxes.add((r, c))
            elif cell == 'X':
                goals.add((r, c))
            elif cell == '$':
                boxes.add((r, c))
                goals.add((r, c))
            elif cell == '%':
                player_pos = (r, c)
                goals.add((r, c))
    
    return grid, player_pos, boxes, goals

# Check if the puzzle is solved
def is_solved(boxes, goals):
    return boxes == goals

# Get possible moves for the player
def get_moves(grid, player_pos, boxes):
    moves = []
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    
    for move, (dr, dc) in directions.items():
        new_r, new_c = player_pos[0] + dr, player_pos[1] + dc
        if grid[new_r][new_c] == '+':
            continue  # Wall
        if (new_r, new_c) in boxes:
            # Check if the box can be pushed
            box_new_r, box_new_c = new_r + dr, new_c + dc
            if grid[box_new_r][box_new_c] == '+' or (box_new_r, box_new_c) in boxes:
                continue  # Wall or another box
            moves.append((move, (new_r, new_c), (box_new_r, box_new_c)))
        else:
            moves.append((move, (new_r, new_c), None))
    
    return moves

# Perform a breadth-first search to find the solution
def solve_sokoban(initial_state):
    grid, player_pos, boxes, goals = parse_state(initial_state)
    queue = deque([(player_pos, boxes, "")])
    visited = set()
    
    while queue:
        player_pos, boxes, path = queue.popleft()
        
        if is_solved(boxes, goals):
            return path
        
        for move, new_player_pos, new_box_pos in get_moves(grid, player_pos, boxes):
            new_boxes = set(boxes)
            if new_box_pos:
                new_boxes.remove((new_player_pos[0], new_player_pos[1]))
                new_boxes.add(new_box_pos)
            
            state = (new_player_pos, frozenset(new_boxes))
            if state not in visited:
                visited.add(state)
                queue.append((new_player_pos, new_boxes, path + move))
    
    return "No solution"

# Find the solution
solution = solve_sokoban(initial_state)
print(solution)