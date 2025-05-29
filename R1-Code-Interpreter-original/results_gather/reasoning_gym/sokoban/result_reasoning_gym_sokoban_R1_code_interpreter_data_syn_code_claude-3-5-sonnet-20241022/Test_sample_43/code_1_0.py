from collections import deque

def is_valid(x, y, width, height, walls):
    return 0 <= x < height and 0 <= y < width and not walls[x][y]

def get_next_pos(x, y, dx, dy):
    return x + dx, y + dy

def serialize_state(player_pos, boxes):
    return (player_pos, tuple(sorted(boxes)))

def is_goal_state(boxes, goals):
    return set(boxes) == set(goals)

def solve_sokoban(walls, player_pos, boxes, goals):
    directions = {'L': (0, -1), 'R': (0, 1), 'U': (-1, 0), 'D': (1, 0)}
    height, width = len(walls), len(walls[0])
    visited = set()
    queue = deque([(player_pos, frozenset(boxes), "")])
    
    while queue:
        curr_player, curr_boxes, path = queue.popleft()
        state = serialize_state(curr_player, curr_boxes)
        
        if state in visited:
            continue
        visited.add(state)
        
        if is_goal_state(curr_boxes, goals):
            return path
        
        for move, (dx, dy) in directions.items():
            new_player_x, new_player_y = get_next_pos(curr_player[0], curr_player[1], dx, dy)
            
            if not is_valid(new_player_x, new_player_y, width, height, walls):
                continue
                
            if (new_player_x, new_player_y) in curr_boxes:
                new_box_x, new_box_y = get_next_pos(new_player_x, new_player_y, dx, dy)
                
                if not is_valid(new_box_x, new_box_y, width, height, walls) or \
                   (new_box_x, new_box_y) in curr_boxes:
                    continue
                    
                new_boxes = set(curr_boxes)
                new_boxes.remove((new_player_x, new_player_y))
                new_boxes.add((new_box_x, new_box_y))
                queue.append(((new_player_x, new_player_y), frozenset(new_boxes), path + move))
            else:
                queue.append(((new_player_x, new_player_y), frozenset(curr_boxes), path + move))
    
    return None

# Initialize the puzzle
walls = [
    [1,1,1,1,1,1,1,1],
    [1,0,0,0,0,0,0,1],
    [1,0,0,0,1,1,0,1],
    [1,0,1,1,0,1,0,1],
    [1,1,1,0,1,1,0,1],
    [1,1,1,1,1,1,1,1]
]

# Initial state
player_pos = (1, 5)  # * position
boxes = {(2,4), (2,5), (4,2), (3,5)}  # @ and $ positions
goals = {(3,2), (3,3), (4,1), (4,3)}  # X positions

solution = solve_sokoban(walls, player_pos, boxes, goals)
print(solution if solution else "No solution found")