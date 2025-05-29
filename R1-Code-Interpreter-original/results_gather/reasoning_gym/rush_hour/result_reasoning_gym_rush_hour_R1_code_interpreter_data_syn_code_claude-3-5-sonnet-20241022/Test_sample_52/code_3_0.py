from heapq import heappush, heappop
from collections import defaultdict

def parse_board(board_str):
    return [list(row) for row in board_str.split('\n') if row]

def find_vehicles(board):
    vehicles = defaultdict(list)
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] not in ['.', 'x']:
                vehicles[board[i][j]].append((i, j))
    return {k: sorted(v) for k, v in vehicles.items()}

def get_orientation(coords):
    return 'horizontal' if coords[0][0] == coords[1][0] else 'vertical'

def manhattan_distance_to_exit(coords, board_width):
    # For red car (A), calculate distance to exit
    if get_orientation(coords) == 'horizontal':
        return board_width - 1 - max(y for _, y in coords)
    return float('inf')

def is_valid_move(board, coords, steps):
    orientation = get_orientation(coords)
    if orientation == 'horizontal':
        row = coords[0][0]
        min_col = min(c[1] for c in coords)
        max_col = max(c[1] for c in coords)
        new_min = min_col + steps
        new_max = max_col + steps
        if new_min < 0 or new_max >= len(board[0]):
            return False
        check_range = range(min_col, new_min) if steps < 0 else range(max_col + 1, new_max + 1)
        return all(board[row][j] == '.' for j in check_range)
    else:
        col = coords[0][1]
        min_row = min(c[0] for c in coords)
        max_row = max(c[0] for c in coords)
        new_min = min_row + steps
        new_max = max_row + steps
        if new_min < 0 or new_max >= len(board):
            return False
        check_range = range(min_row, new_min) if steps < 0 else range(max_row + 1, new_max + 1)
        return all(board[i][col] == '.' for i in check_range)

def make_move(board, vehicle, coords, steps):
    new_board = [row[:] for row in board]
    orientation = get_orientation(coords)
    # Clear old positions
    for x, y in coords:
        new_board[x][y] = '.'
    # Set new positions
    for x, y in coords:
        new_x = x + (steps if orientation == 'vertical' else 0)
        new_y = y + (steps if orientation == 'horizontal' else 0)
        new_board[new_x][new_y] = vehicle
    return new_board

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_puzzle():
    initial_board = parse_board("BBH...\nx.HCCC\n.AAI.K\n.G.IJK\n.GDDJK\n.EEEFF")
    vehicles = find_vehicles(initial_board)
    
    # Priority queue: (priority, moves_count, board_state, moves_list)
    queue = [(0, 0, initial_board, [])]
    seen = {board_to_string(initial_board)}
    
    while queue:
        _, moves_count, current_board, moves = heappop(queue)
        vehicles = find_vehicles(current_board)
        
        # Check if solved
        red_car = vehicles['A']
        if max(y for _, y in red_car) == len(current_board[0]) - 1:
            return moves
        
        # Try moves for each vehicle
        for vehicle, coords in vehicles.items():
            # Calculate maximum possible movement range
            max_range = 3 if vehicle == 'A' else 2
            for steps in range(-max_range, max_range + 1):
                if steps == 0:
                    continue
                    
                if is_valid_move(current_board, coords, steps):
                    new_board = make_move(current_board, vehicle, coords, steps)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in seen:
                        seen.add(board_str)
                        new_moves = moves + [f"{vehicle}{'+' if steps > 0 else ''}{steps}"]
                        
                        # Calculate priority based on red car's position and moves count
                        new_vehicles = find_vehicles(new_board)
                        h_score = manhattan_distance_to_exit(new_vehicles['A'], len(current_board[0]))
                        priority = h_score + len(new_moves)
                        
                        heappush(queue, (priority, len(new_moves), new_board, new_moves))
    
    return None

# Solve and print the solution
solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")