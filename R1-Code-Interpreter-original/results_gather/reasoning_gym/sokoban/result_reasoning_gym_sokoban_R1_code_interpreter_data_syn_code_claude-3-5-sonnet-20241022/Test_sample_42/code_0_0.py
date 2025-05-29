from collections import deque
import copy

def get_initial_state():
    puzzle = [
        "++++++++++"
        "+--X---@*+"
        "+---@X@@@+"
        "+---X+X$-+"
        "++-++$++X+"
        "+X--++++$+"
        "+-@-+++++"
        "+----++++"
        "+----+++$+"
        "++++++++++"
    ]
    return ''.join(puzzle)

def is_goal_state(state):
    # Count boxes not on goals
    boxes = state.count('@')
    return boxes == 0

def get_player_pos(state):
    width = 10
    if '*' in state:
        pos = state.index('*')
    else:
        pos = state.index('@')
    return (pos // width, pos % width)

def get_valid_moves(state):
    width = 10
    height = 10
    grid = [list(state[i:i+width]) for i in range(0, len(state), width)]
    player_row, player_col = get_player_pos(state)
    moves = []
    
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    
    for direction, (dr, dc) in directions.items():
        new_row, new_col = player_row + dr, player_col + dc
        
        if 0 <= new_row < height and 0 <= new_col < width:
            curr_cell = grid[player_row][player_col]
            next_cell = grid[new_row][new_col]
            
            # Simple move
            if next_cell in '-X':
                new_grid = copy.deepcopy(grid)
                new_grid[player_row][player_col] = '-' if curr_cell == '@' else 'X'
                new_grid[new_row][new_col] = '@' if next_cell == '-' else '*'
                moves.append((direction, ''.join(''.join(row) for row in new_grid)))
            
            # Push box
            elif next_cell in '@$':
                push_row, push_col = new_row + dr, new_col + dc
                if 0 <= push_row < height and 0 <= push_col < width:
                    push_cell = grid[push_row][push_col]
                    if push_cell in '-X':
                        new_grid = copy.deepcopy(grid)
                        new_grid[player_row][player_col] = '-' if curr_cell == '@' else 'X'
                        new_grid[new_row][new_col] = '@' if next_cell == '@' else '*'
                        new_grid[push_row][push_col] = '@' if push_cell == '-' else '$'
                        moves.append((direction, ''.join(''.join(row) for row in new_grid)))
    
    return moves

def solve_sokoban():
    initial_state = get_initial_state()
    queue = deque([(initial_state, "")])
    visited = {initial_state}
    
    while queue:
        state, path = queue.popleft()
        
        if is_goal_state(state):
            return path
        
        for move, new_state in get_valid_moves(state):
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + move))
    
    return None

solution = solve_sokoban()
print(solution if solution else "No solution found")