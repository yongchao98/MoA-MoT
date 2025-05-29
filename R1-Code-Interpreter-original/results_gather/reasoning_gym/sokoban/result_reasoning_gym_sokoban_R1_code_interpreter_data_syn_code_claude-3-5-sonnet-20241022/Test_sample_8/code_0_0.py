from collections import deque
import copy

# Initialize the puzzle
puzzle = [
    "+++++++++",
    "+*$---+++",
    "+@----+++",
    "+X@@---++",
    "++X$-X-++",
    "+++X-@@++",
    "+$+X----+",
    "++-@----+",
    "++-----X+",
    "+++++++++"
]

def get_initial_state():
    player_pos = None
    boxes = set()
    goals = set()
    
    for i in range(len(puzzle)):
        for j in range(len(puzzle[i])):
            if puzzle[i][j] in ['*', '@']:
                player_pos = (i, j)
            if puzzle[i][j] in ['@', '$']:
                boxes.add((i, j))
            if puzzle[i][j] in ['X', '$', '%']:
                goals.add((i, j))
    return player_pos, boxes, goals

def is_valid_move(pos, boxes):
    return (0 <= pos[0] < len(puzzle) and 
            0 <= pos[1] < len(puzzle[0]) and 
            puzzle[pos[0]][pos[1]] != '+' and 
            pos not in boxes)

def get_moves():
    return {
        'U': (-1, 0, 'U'),
        'D': (1, 0, 'D'),
        'L': (0, -1, 'L'),
        'R': (0, 1, 'R')
    }

def solve_sokoban():
    initial_player, initial_boxes, goals = get_initial_state()
    queue = deque([(initial_player, initial_boxes, "")])
    visited = set()
    
    while queue:
        player, boxes, path = queue.popleft()
        state = (player, tuple(sorted(boxes)))
        
        if state in visited:
            continue
        visited.add(state)
        
        if boxes == goals:
            return path
            
        for move, (dy, dx, direction) in get_moves().items():
            new_player = (player[0] + dy, player[1] + dx)
            
            if not is_valid_move(new_player, boxes):
                continue
                
            new_boxes = set(boxes)
            if new_player in boxes:
                box_new_pos = (new_player[0] + dy, new_player[1] + dx)
                if not is_valid_move(box_new_pos, boxes - {new_player}):
                    continue
                new_boxes.remove(new_player)
                new_boxes.add(box_new_pos)
            
            queue.append((new_player, new_boxes, path + direction))
    
    return None

solution = solve_sokoban()
print(solution if solution else "No solution found")