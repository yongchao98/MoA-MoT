from heapq import heappush, heappop

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def is_deadlock(grid, box_pos, goals):
    # Simple deadlock detection for corners
    if box_pos not in goals:
        row, col = box_pos
        if (grid[row-1][col] == '+' and grid[row][col-1] == '+') or \
           (grid[row-1][col] == '+' and grid[row][col+1] == '+') or \
           (grid[row+1][col] == '+' and grid[row][col-1] == '+') or \
           (grid[row+1][col] == '+' and grid[row][col+1] == '+'):
            return True
    return False

def solve_sokoban(puzzle):
    # Convert puzzle to grid
    grid = [list(row.strip()) for row in puzzle]
    rows, cols = len(grid), len(grid[0])
    
    # Find initial state
    goals = set()
    boxes = set()
    player = None
    
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] == 'X':
                goals.add((i, j))
            elif grid[i][j] == '$':
                goals.add((i, j))
                boxes.add((i, j))
            elif grid[i][j] == '@':
                boxes.add((i, j))
            elif grid[i][j] == '*':
                player = (i, j)
    
    # Directions: Right, Left, Down, Up
    moves = [(0, 1, 'R'), (0, -1, 'L'), (1, 0, 'D'), (-1, 0, 'U')]
    visited = set()
    queue = [(0, player, frozenset(boxes), "")]
    visited.add((player, frozenset(boxes)))
    
    max_steps = 100  # Limit search depth
    
    while queue and len(queue[0][3]) < max_steps:
        _, current_player, current_boxes, path = heappop(queue)
        
        if all(box in goals for box in current_boxes):
            return path
        
        for dx, dy, move in moves:
            new_player = (current_player[0] + dx, current_player[1] + dy)
            
            # Check if move is valid
            if grid[new_player[0]][new_player[1]] == '+':
                continue
                
            new_boxes = set(current_boxes)
            if new_player in current_boxes:
                new_box_pos = (new_player[0] + dx, new_player[1] + dy)
                
                # Check if box can be pushed
                if grid[new_box_pos[0]][new_box_pos[1]] == '+' or \
                   new_box_pos in new_boxes or \
                   is_deadlock(grid, new_box_pos, goals):
                    continue
                    
                new_boxes.remove(new_player)
                new_boxes.add(new_box_pos)
            
            new_boxes_frozen = frozenset(new_boxes)
            if (new_player, new_boxes_frozen) not in visited:
                visited.add((new_player, new_boxes_frozen))
                # Simple heuristic: distance to closest goal for each box
                h = sum(min(manhattan_distance(box, goal) for goal in goals) for box in new_boxes)
                heappush(queue, (len(path) + h, new_player, new_boxes_frozen, path + move))
    
    return None

puzzle = [
    "++++++++++",
    "+X-@-$*+++",
    "+-@X--@-++",
    "+-X-@---++",
    "+-X--@-X++",
    "+-@---++++",
    "+-@X-@---+",
    "++X-$----+",
    "+++++X---+",
    "++++++++++"
]

solution = solve_sokoban(puzzle)
print(solution)