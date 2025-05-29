from collections import deque
import copy

def is_valid(grid, row, col):
    return 0 <= row < len(grid) and 0 <= col < len(grid[0]) and grid[row][col] != '+'

def get_next_positions(grid, player_pos):
    moves = {'U': (-1, 0, 'U'), 'D': (1, 0, 'D'), 'L': (0, -1, 'L'), 'R': (0, 1, 'R')}
    valid_moves = []
    
    for move, (dr, dc, direction) in moves.items():
        new_r, new_c = player_pos[0] + dr, player_pos[1] + dc
        
        if is_valid(grid, new_r, new_c):
            if grid[new_r][new_c] in '-X':
                valid_moves.append((new_r, new_c, direction, False))
            elif grid[new_r][new_c] in '@$':
                box_r, box_c = new_r + dr, new_c + dc
                if is_valid(grid, box_r, box_c) and grid[box_r][box_c] in '-X':
                    valid_moves.append((new_r, new_c, direction, True))
    
    return valid_moves

def is_goal_state(grid):
    for row in grid:
        for cell in row:
            if cell == '@':
                return False
    return True

def state_to_string(grid):
    return ''.join(''.join(row) for row in grid)

def solve_sokoban(initial_grid):
    rows, cols = len(initial_grid), len(initial_grid[0])
    
    # Find player position
    player_pos = None
    for i in range(rows):
        for j in range(cols):
            if initial_grid[i][j] in '*%':
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    queue = deque([(initial_grid, player_pos, "")])
    visited = {state_to_string(initial_grid)}
    
    while queue:
        current_grid, current_pos, path = queue.popleft()
        
        if is_goal_state(current_grid):
            return path
        
        for new_r, new_c, direction, is_push in get_next_positions(current_grid, current_pos):
            new_grid = [list(row) for row in current_grid]
            r, c = current_pos
            
            # Update player position
            if new_grid[r][c] == '%':
                new_grid[r][c] = 'X'
            else:
                new_grid[r][c] = '-'
            
            if new_grid[new_r][new_c] == 'X':
                new_grid[new_r][new_c] = '%'
            else:
                new_grid[new_r][new_c] = '*'
            
            # If pushing a box
            if is_push:
                dr = new_r - r
                dc = new_c - c
                box_r, box_c = new_r + dr, new_c + dc
                
                if new_grid[box_r][box_c] == 'X':
                    new_grid[box_r][box_c] = '$'
                else:
                    new_grid[box_r][box_c] = '@'
            
            new_state = state_to_string(new_grid)
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_grid, (new_r, new_c), path + direction))

    return None

# Initialize the puzzle
initial_grid = [
    list("++++++"),
    list("+-@X+"),
    list("+%@-+"),
    list("+-$-+"),
    list("+---+"),
    list("++X@+"),
    list("++++++")
]

solution = solve_sokoban(initial_grid)
print(solution if solution else "No solution found")