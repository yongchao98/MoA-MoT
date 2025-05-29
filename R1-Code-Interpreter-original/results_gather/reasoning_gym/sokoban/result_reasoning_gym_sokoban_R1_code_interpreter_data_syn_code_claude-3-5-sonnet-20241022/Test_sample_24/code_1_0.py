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
    
    def dfs(player, boxes, path, depth, visited):
        if depth > 30:  # Strict depth limit
            return None
        
        if all(box in goals for box in boxes):
            return path
        
        # Directions: prioritize moving towards goals
        moves = [(0, 1, 'R'), (0, -1, 'L'), (1, 0, 'D'), (-1, 0, 'U')]
        
        for dx, dy, move in moves:
            new_player = (player[0] + dx, player[1] + dy)
            
            # Basic validity checks
            if grid[new_player[0]][new_player[1]] == '+':
                continue
                
            new_boxes = set(boxes)
            if new_player in boxes:
                new_box_pos = (new_player[0] + dx, new_player[1] + dy)
                
                # Check if box can be pushed
                if (grid[new_box_pos[0]][new_box_pos[1]] == '+' or 
                    new_box_pos in new_boxes):
                    continue
                    
                new_boxes.remove(new_player)
                new_boxes.add(new_box_pos)
            
            state = (new_player, tuple(sorted(new_boxes)))
            if state not in visited:
                visited.add(state)
                result = dfs(new_player, new_boxes, path + move, depth + 1, visited)
                if result:
                    return result
        
        return None

    # Start DFS with initial state
    initial_state = (player, tuple(sorted(boxes)))
    visited = {initial_state}
    solution = dfs(player, boxes, "", 0, visited)
    
    return solution

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