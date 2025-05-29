from collections import deque

def is_valid(grid, row, col):
    return 0 <= row < len(grid) and 0 <= col < len(grid[0]) and grid[row][col] != '+'

def get_initial_state(grid):
    players = []
    boxes = set()
    goals = set()
    
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] in ['@', '*', '%']:  # All possible player positions
                players.append((i, j))
            if grid[i][j] in ['@', '$']:  # All boxes
                boxes.add((i, j))
            if grid[i][j] in ['X', '$', '%']:  # All goals
                goals.add((i, j))
    
    return players, boxes, goals

def get_state_key(players, boxes):
    return (tuple(sorted(players)), tuple(sorted(boxes)))

def solve_sokoban(grid):
    players, boxes, goals = get_initial_state(grid)
    if not players or not boxes or not goals:
        return None
    
    directions = [('U', -1, 0), ('D', 1, 0), ('L', 0, -1), ('R', 0, 1)]
    initial_state = (players, boxes, "")
    queue = deque([initial_state])
    visited = set()
    
    while queue:
        curr_players, curr_boxes, path = queue.popleft()
        
        if curr_boxes == goals:
            return path
        
        state_key = get_state_key(curr_players, curr_boxes)
        if state_key in visited:
            continue
        visited.add(state_key)
        
        # Try moving each player
        for player in curr_players:
            for move, dy, dx in directions:
                new_player_pos = (player[0] + dy, player[1] + dx)
                
                if not is_valid(grid, new_player_pos[0], new_player_pos[1]):
                    continue
                
                new_boxes = set(curr_boxes)
                new_players = [p for p in curr_players if p != player]
                
                if new_player_pos in curr_boxes:  # If we're pushing a box
                    new_box_pos = (new_player_pos[0] + dy, new_player_pos[1] + dx)
                    
                    if not is_valid(grid, new_box_pos[0], new_box_pos[1]) or \
                       new_box_pos in new_boxes:
                        continue
                    
                    new_boxes.remove(new_player_pos)
                    new_boxes.add(new_box_pos)
                
                new_players.append(new_player_pos)
                queue.append((new_players, new_boxes, path + move))
    
    return None

# Initialize the puzzle
puzzle = [
    list("+++++++" ),
    list("++$++++" ),
    list("+X+++$+" ),
    list("+-+++++" ),
    list("+---+$+" ),
    list("+--@+$+" ),
    list("+--X+++" ),
    list("+@--+++" ),
    list("+*@-X++" ),
    list("+++++++" )
]

solution = solve_sokoban(puzzle)
print(solution if solution else "No solution found")