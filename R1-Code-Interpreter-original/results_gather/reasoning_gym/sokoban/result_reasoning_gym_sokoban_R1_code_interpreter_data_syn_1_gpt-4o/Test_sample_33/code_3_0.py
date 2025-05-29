from heapq import heappop, heappush

# Define the initial state of the puzzle
initial_state = [
    "+ + + + + + + + +",
    "+ - - - - - - - +",
    "+ X - @ X - @ - +",
    "+ @ X @ - @ @ X +",
    "+ - - - - X - - +",
    "+ - * @ X - - - +",
    "+ - - - X - - - +",
    "+ + + + + + + + +"
]

# Convert the initial state into a more manageable format
def parse_state(state):
    grid = [list(row.replace(" ", "")) for row in state]
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

# Check if the current state is a winning state
def is_winning(boxes, goals):
    return boxes == goals

# Get possible moves from the current position
def get_moves(grid, player_pos, boxes):
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    moves = []
    
    for move, (dr, dc) in directions.items():
        new_r, new_c = player_pos[0] + dr, player_pos[1] + dc
        if grid[new_r][new_c] in '-X' and (new_r, new_c) not in boxes:
            moves.append((move, (new_r, new_c), boxes))
        elif (new_r, new_c) in boxes:
            box_new_r, box_new_c = new_r + dr, new_c + dc
            if grid[box_new_r][box_new_c] in '-X' and (box_new_r, box_new_c) not in boxes:
                new_boxes = boxes - {(new_r, new_c)} | {(box_new_r, box_new_c)}
                moves.append((move, (new_r, new_c), new_boxes))
    
    return moves

# Heuristic function for A* search
def heuristic(boxes, goals):
    return sum(min(abs(bx - gx) + abs(by - gy) for gx, gy in goals) for bx, by in boxes)

# Perform A* search to find the solution
def solve_sokoban(initial_state):
    grid, player_pos, boxes, goals = parse_state(initial_state)
    queue = [(0 + heuristic(boxes, goals), 0, player_pos, boxes, "")]
    visited = set()
    
    while queue:
        _, cost, player_pos, boxes, path = heappop(queue)
        
        if is_winning(boxes, goals):
            return path
        
        for move, new_player_pos, new_boxes in get_moves(grid, player_pos, boxes):
            state = (new_player_pos, frozenset(new_boxes))
            if state not in visited:
                visited.add(state)
                heappush(queue, (cost + 1 + heuristic(new_boxes, goals), cost + 1, new_player_pos, new_boxes, path + move))
    
    return "No solution"

# Solve the puzzle
solution = solve_sokoban(initial_state)
print(solution)