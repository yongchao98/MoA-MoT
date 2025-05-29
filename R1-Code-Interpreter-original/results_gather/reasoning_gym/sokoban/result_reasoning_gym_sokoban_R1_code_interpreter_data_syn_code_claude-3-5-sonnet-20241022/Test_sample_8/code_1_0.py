from collections import deque

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
            if puzzle[i][j] == '*':  # player
                player_pos = (i, j)
            elif puzzle[i][j] == '@':  # box
                boxes.add((i, j))
            elif puzzle[i][j] == 'X':  # goal
                goals.add((i, j))
            elif puzzle[i][j] == '$':  # box on goal
                boxes.add((i, j))
                goals.add((i, j))
            elif puzzle[i][j] == '%':  # player on goal
                player_pos = (i, j)
                goals.add((i, j))
    return player_pos, boxes, goals

def is_valid_pos(pos):
    return (0 <= pos[0] < len(puzzle) and 
            0 <= pos[1] < len(puzzle[0]) and 
            puzzle[pos[0]][pos[1]] != '+')

def solve_sokoban():
    initial_player, initial_boxes, goals = get_initial_state()
    queue = deque([(initial_player, frozenset(initial_boxes), "")])
    visited = set()
    moves = [('U', -1, 0), ('D', 1, 0), ('L', 0, -1), ('R', 0, 1)]
    
    while queue:
        player, boxes, path = queue.popleft()
        
        if boxes == goals:
            return path
        
        state = (player, boxes)
        if state in visited:
            continue
        visited.add(state)
        
        for direction, dy, dx in moves:
            new_player = (player[0] + dy, player[1] + dx)
            
            if not is_valid_pos(new_player):
                continue
                
            if new_player in boxes:
                # Try to push the box
                new_box_pos = (new_player[0] + dy, new_player[1] + dx)
                if not is_valid_pos(new_box_pos) or new_box_pos in boxes:
                    continue
                    
                # Create new box configuration
                new_boxes = set(boxes)
                new_boxes.remove(new_player)
                new_boxes.add(new_box_pos)
                queue.append((new_player, frozenset(new_boxes), path + direction))
            else:
                # Simple player move
                queue.append((new_player, boxes, path + direction))
    
    return None

solution = solve_sokoban()
print(solution if solution else "No solution found")