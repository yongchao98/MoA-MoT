def get_player_pos(grid):
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            if grid[i][j] in ['*', '%']:
                return (i, j)
    return None

def is_valid_move(grid, y, x):
    return 0 <= y < len(grid) and 0 <= x < len(grid[0]) and grid[y][x] != '+'

def try_move(grid, player_pos, dy, dx, moves_limit=25):
    moves = []
    visited = set()
    stack = [(player_pos, grid, [])]
    
    while stack and len(moves) < moves_limit:
        pos, current_grid, path = stack.pop()
        state_key = str(current_grid)
        
        if state_key in visited:
            continue
        visited.add(state_key)
        
        y, x = pos
        new_y, new_x = y + dy, x + dx
        
        if not is_valid_move(current_grid, new_y, new_x):
            continue
            
        new_grid = [row[:] for row in current_grid]
        
        # Moving to empty space or goal
        if current_grid[new_y][new_x] in '-X':
            new_grid[y][x] = 'X' if current_grid[y][x] == '*' else '-'
            new_grid[new_y][new_x] = '%' if current_grid[new_y][new_x] == 'X' else '*'
            moves.append(path + ['UDLR'[['U', 'D', 'L', 'R'].index((dy, dx))]])
            
        # Moving box
        elif current_grid[new_y][new_x] in '@$':
            push_y, push_x = new_y + dy, new_x + dx
            if is_valid_move(current_grid, push_y, push_x) and current_grid[push_y][push_x] in '-X':
                new_grid[y][x] = 'X' if current_grid[y][x] == '*' else '-'
                new_grid[new_y][new_x] = '%' if current_grid[new_y][new_x] == '$' else '*'
                new_grid[push_y][push_x] = '$' if current_grid[push_y][push_x] == 'X' else '@'
                moves.append(path + ['UDLR'[['U', 'D', 'L', 'R'].index((dy, dx))]])
                
    return moves[0] if moves else None

def solve_simple(grid):
    solution = []
    max_moves = 25
    directions = [(-1, 0, 'U'), (1, 0, 'D'), (0, -1, 'L'), (0, 1, 'R')]
    
    while len(solution) < max_moves:
        player_pos = get_player_pos(grid)
        if not player_pos:
            break
            
        for dy, dx, dir_char in directions:
            move = try_move(grid, player_pos, dy, dx)
            if move:
                solution.extend(move)
                break
                
        if not move:
            break
            
    return ''.join(solution) if solution else None

# Initialize puzzle
puzzle = [
    list("+++++++++++"),
    list("+---------+"),
    list("+-@XX--@--+"),
    list("+-$-@X---+"),
    list("+--X-----+"),
    list("+-$@+--+-+"),
    list("+---X+---+"),
    list("+$@@X--@-+"),
    list("+*----+X-+"),
    list("+++++++++++")
]

solution = solve_simple(puzzle)
print(solution)