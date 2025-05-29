from collections import defaultdict, deque

def parse_board(board_str):
    board = [list(row) for row in board_str.split('\n') if row]
    return board

def find_vehicles(board):
    vehicles = defaultdict(list)
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] not in ['.', 'x']:
                vehicles[board[i][j]].append((i, j))
    for v in vehicles:
        vehicles[v].sort()
    return vehicles

def get_orientation(coords):
    return 'horizontal' if coords[0][0] == coords[1][0] else 'vertical'

def is_valid_move(board, vehicle, coords, steps):
    orientation = get_orientation(coords)
    if orientation == 'horizontal':
        row = coords[0][0]
        min_col = min(c[1] for c in coords)
        max_col = max(c[1] for c in coords)
        if steps > 0 and max_col + steps >= len(board[0]): return False
        if steps < 0 and min_col + steps < 0: return False
        check_range = range(max_col + 1, max_col + steps + 1) if steps > 0 else range(min_col + steps, min_col)
        return all(board[row][j] in ['.', vehicle] for j in check_range)
    else:
        col = coords[0][1]
        min_row = min(c[0] for c in coords)
        max_row = max(c[0] for c in coords)
        if steps > 0 and max_row + steps >= len(board): return False
        if steps < 0 and min_row + steps < 0: return False
        check_range = range(max_row + 1, max_row + steps + 1) if steps > 0 else range(min_row + steps, min_row)
        return all(board[i][col] in ['.', vehicle] for i in check_range)

def make_move(board, vehicle, coords, steps):
    new_board = [row[:] for row in board]
    orientation = get_orientation(coords)
    for x, y in coords:
        new_board[x][y] = '.'
    for x, y in coords:
        new_x = x + (steps if orientation == 'vertical' else 0)
        new_y = y + (steps if orientation == 'horizontal' else 0)
        new_board[new_x][new_y] = vehicle
    return new_board

def solve_puzzle():
    initial_board = parse_board("BBH...\nx.HCCC\n.AAI.K\n.G.IJK\n.GDDJK\n.EEEFF")
    queue = deque([(initial_board, [])])
    seen = {str(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        vehicles = find_vehicles(current_board)
        
        if max(y for _, y in vehicles['A']) == len(current_board[0]) - 1:
            return moves
            
        # Prioritize moves that might lead to solution
        for vehicle in sorted(vehicles.keys(), key=lambda x: (x != 'A', x)):
            coords = vehicles[vehicle]
            orientation = get_orientation(coords)
            
            # Try larger moves first for potential shortcuts
            steps_range = range(-3, 4) if vehicle == 'A' else range(-2, 3)
            for steps in steps_range:
                if steps == 0: continue
                if is_valid_move(current_board, vehicle, coords, steps):
                    new_board = make_move(current_board, vehicle, coords, steps)
                    board_str = str(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        new_moves = moves + [f"{vehicle}{'+' if steps > 0 else ''}{steps}"]
                        queue.append((new_board, new_moves))
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")