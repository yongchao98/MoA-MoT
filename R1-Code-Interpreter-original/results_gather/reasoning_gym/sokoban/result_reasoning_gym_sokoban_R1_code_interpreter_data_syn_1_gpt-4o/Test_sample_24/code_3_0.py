from heapq import heappop, heappush

# Define the initial state of the puzzle
initial_state = [
    "+ + + + + + + + + +",
    "+ X - @ - $ * + + +",
    "+ - @ X - - @ - + +",
    "+ - X - @ - - - + +",
    "+ - X - - @ - X + +",
    "+ - @ - - - - + + +",
    "+ - @ X - @ - - - +",
    "+ + X - $ - - - - +",
    "+ + + + + X - - - +",
    "+ + + + + + + + + +"
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
            elif cell == '$':
                boxes.add((r, c))
                goals.add((r, c))
            elif cell == 'X':
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
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    moves = []
    
    for move, (dr, dc) in directions.items():
        new_r, new_c = player_pos[0] + dr, player_pos[1] + dc
        if grid[new_r][new_c] in ('-', 'X'):
            moves.append((move, (new_r, new_c)))
        elif grid[new_r][new_c] in ('@', '$'):
            box_r, box_c = new_r + dr, new_c + dc
            if (box_r, box_c) not in boxes and grid[box_r][box_c] in ('-', 'X'):
                moves.append((move, (new_r, new_c)))
    
    return moves

# Heuristic function to estimate the cost to reach the goal
def heuristic(boxes, goals):
    return sum(min(abs(bx - gx) + abs(by - gy) for gx, gy in goals) for bx, by in boxes)

# Check for deadlocks
def is_deadlock(boxes, grid):
    for box in boxes:
        r, c = box
        if grid[r][c] == 'X':
            continue
        # Check boundaries before accessing grid
        if (r > 0 and r < len(grid) - 1 and c > 0 and c < len(grid[0]) - 1):
            if (grid[r-1][c] == '+' or grid[r+1][c] == '+') and (grid[r][c-1] == '+' or grid[r][c+1] == '+'):
                return True
    return False

# Perform an A* search to find the solution
def solve_sokoban(initial_state):
    grid, player_pos, boxes, goals = parse_state(initial_state)
    queue = [(0 + heuristic(boxes, goals), 0, player_pos, boxes, "")]
    visited = set()
    
    while queue:
        _, cost, player_pos, boxes, path = heappop(queue)
        
        if is_solved(boxes, goals):
            return path
        
        if is_deadlock(boxes, grid):
            continue
        
        for move, new_player_pos in get_moves(grid, player_pos, boxes):
            new_boxes = set(boxes)
            if new_player_pos in boxes:
                dr, dc = new_player_pos[0] - player_pos[0], new_player_pos[1] - player_pos[1]
                new_box_pos = (new_player_pos[0] + dr, new_player_pos[1] + dc)
                new_boxes.remove(new_player_pos)
                new_boxes.add(new_box_pos)
            
            state = (new_player_pos, frozenset(new_boxes))
            if state not in visited:
                visited.add(state)
                new_cost = cost + 1
                priority = new_cost + heuristic(new_boxes, goals)
                heappush(queue, (priority, new_cost, new_player_pos, new_boxes, path + move))
    
    return "No solution"

# Find the solution to the puzzle
solution = solve_sokoban(initial_state)
print(solution)