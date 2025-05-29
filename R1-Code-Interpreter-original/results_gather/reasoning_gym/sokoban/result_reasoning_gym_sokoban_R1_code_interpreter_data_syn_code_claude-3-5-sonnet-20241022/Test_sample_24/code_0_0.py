from heapq import heappush, heappop
from copy import deepcopy

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def get_state_string(boxes, player):
    return str(sorted(boxes)) + str(player)

def is_valid_move(grid, player, boxes, new_player_pos):
    if grid[new_player_pos[0]][new_player_pos[1]] == '+':
        return False
    
    if new_player_pos in boxes:
        box_pos = (new_player_pos[0] + (new_player_pos[0] - player[0]),
                  new_player_pos[1] + (new_player_pos[1] - player[1]))
        if box_pos in boxes or grid[box_pos[0]][box_pos[1]] == '+':
            return False
    return True

def solve_sokoban(grid):
    rows, cols = len(grid), len(grid[0])
    goals = set()
    boxes = set()
    player = None
    
    # Parse initial state
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] in ['X', '$', '%']:
                goals.add((i, j))
            if grid[i][j] in ['@', '$']:
                boxes.add((i, j))
            if grid[i][j] in ['*', '%']:
                player = (i, j)

    moves = [(0, 1, 'R'), (0, -1, 'L'), (1, 0, 'D'), (-1, 0, 'U')]
    visited = set()
    queue = [(0, 0, boxes, player, "")]
    visited.add(get_state_string(boxes, player))
    
    while queue:
        _, cost, current_boxes, current_player, path = heappop(queue)
        
        if all(box in goals for box in current_boxes):
            return path
            
        for dx, dy, move in moves:
            new_player = (current_player[0] + dx, current_player[1] + dy)
            new_boxes = set(current_boxes)
            
            if is_valid_move(grid, current_player, current_boxes, new_player):
                if new_player in current_boxes:
                    new_box_pos = (new_player[0] + dx, new_player[1] + dy)
                    new_boxes.remove(new_player)
                    new_boxes.add(new_box_pos)
                
                state_string = get_state_string(new_boxes, new_player)
                if state_string not in visited:
                    visited.add(state_string)
                    
                    # Calculate heuristic
                    h = sum(min(manhattan_distance(box, goal) for goal in goals) for box in new_boxes)
                    heappush(queue, (cost + 1 + h, cost + 1, new_boxes, new_player, path + move))

# Initialize the puzzle
puzzle = [
    "+++++++++++",
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