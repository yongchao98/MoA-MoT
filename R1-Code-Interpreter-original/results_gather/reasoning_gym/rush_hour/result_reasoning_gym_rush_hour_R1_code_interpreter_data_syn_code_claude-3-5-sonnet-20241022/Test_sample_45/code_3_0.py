from collections import defaultdict

def parse_board(board):
    height, width = len(board), len(board[0])
    vehicles = defaultdict(list)
    grid = [list(row) for row in board]
    
    # Find all vehicles and their positions
    for i in range(height):
        for j in range(width):
            if grid[i][j] not in '.x':
                vehicles[grid[i][j]].append((i, j))
    return grid, vehicles

def get_blocking_chain(grid, vehicles):
    # Get the chain of vehicles blocking the red car
    red_car = vehicles['A']
    red_row = red_car[0][0]
    red_right = max(pos[1] for pos in red_car)
    blocking = set()
    to_check = [(red_row, j) for j in range(red_right + 1, len(grid[0]))]
    
    while to_check:
        i, j = to_check.pop(0)
        if j >= len(grid[0]):
            continue
        if grid[i][j] not in '.x':
            vehicle = grid[i][j]
            if vehicle not in blocking:
                blocking.add(vehicle)
                # Add positions that block this vehicle
                v_pos = vehicles[vehicle]
                if v_pos[0][0] == v_pos[1][0]:  # horizontal
                    row = v_pos[0][0]
                    for col in range(min(p[1] for p in v_pos) - 1, max(p[1] for p in v_pos) + 2):
                        if 0 <= col < len(grid[0]) and grid[row][col] not in '.x':
                            to_check.append((row, col))
                else:  # vertical
                    col = v_pos[0][1]
                    for row in range(min(p[0] for p in v_pos) - 1, max(p[0] for p in v_pos) + 2):
                        if 0 <= row < len(grid) and grid[row][col] not in '.x':
                            to_check.append((row, col))
    return blocking

def get_moves(grid, vehicles):
    height, width = len(grid), len(grid[0])
    moves = []
    
    for vehicle, positions in vehicles.items():
        is_horizontal = positions[0][0] == positions[1][0]
        
        if is_horizontal:
            row = positions[0][0]
            min_col = min(p[1] for p in positions)
            max_col = max(p[1] for p in positions)
            
            # Try move left
            if min_col > 0 and grid[row][min_col - 1] == '.':
                moves.append((vehicle, -1))
            # Try move right
            if max_col < width - 1 and grid[row][max_col + 1] == '.':
                moves.append((vehicle, 1))
        else:
            col = positions[0][1]
            min_row = min(p[0] for p in positions)
            max_row = max(p[0] for p in positions)
            
            # Try move up
            if min_row > 0 and grid[min_row - 1][col] == '.':
                moves.append((vehicle, -1))
            # Try move down
            if max_row < height - 1 and grid[max_row + 1][col] == '.':
                moves.append((vehicle, 1))
    
    return moves

def apply_move(grid, vehicles, move):
    vehicle, direction = move
    new_grid = [row[:] for row in grid]
    positions = vehicles[vehicle]
    is_horizontal = positions[0][0] == positions[1][0]
    
    # Clear current positions
    for i, j in positions:
        new_grid[i][j] = '.'
    
    # Add new positions
    new_positions = []
    for i, j in positions:
        new_i = i + (0 if is_horizontal else direction)
        new_j = j + (direction if is_horizontal else 0)
        new_grid[new_i][new_j] = vehicle
        new_positions.append((new_i, new_j))
    
    return new_grid, new_positions

def solve_puzzle(board, max_depth=50):
    def dfs(grid, vehicles, depth, path, seen):
        if depth == 0:
            return None
        
        # Check if solved
        red_car = vehicles['A']
        if max(pos[1] for pos in red_car) == len(grid[0]) - 1:
            return path
        
        # Get blocking chain and prioritize those moves
        blocking = get_blocking_chain(grid, vehicles)
        moves = get_moves(grid, vehicles)
        moves.sort(key=lambda m: m[0] in blocking, reverse=True)
        
        grid_tuple = tuple(tuple(row) for row in grid)
        if grid_tuple in seen:
            return None
        seen.add(grid_tuple)
        
        for move in moves:
            new_grid, new_positions = apply_move(grid, vehicles, move)
            new_vehicles = vehicles.copy()
            new_vehicles[move[0]] = new_positions
            
            move_str = f"{move[0]}{'+' if move[1] > 0 else ''}{move[1]}"
            result = dfs(new_grid, new_vehicles, depth - 1, path + [move_str], seen)
            if result:
                return result
        
        return None
    
    # Initial state
    grid, vehicles = parse_board(board)
    seen = set()
    
    # Try increasing depths until solution is found
    for depth in range(1, max_depth + 1):
        result = dfs(grid, vehicles, depth, [], seen)
        if result:
            return result
    
    return None

# Initial board
board = [
    "BBH.CC",
    ".GH..K",
    ".GAAJK",
    "DD.IJL",
    "...IxL",
    ".EEFFL"
]

solution = solve_puzzle(board)
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")