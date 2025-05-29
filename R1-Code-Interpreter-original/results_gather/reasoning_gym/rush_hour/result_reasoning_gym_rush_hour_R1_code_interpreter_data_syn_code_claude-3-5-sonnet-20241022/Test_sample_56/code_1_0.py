def get_car_data(board):
    cars = {}
    for i, row in enumerate(board):
        for j, cell in enumerate(row):
            if cell != '.' and cell != 'x':
                if cell not in cars:
                    cars[cell] = {'pos': [], 'orientation': None}
                cars[cell]['pos'].append((i, j))
    
    # Set orientation and sort positions
    for car in cars:
        pos = cars[car]['pos']
        cars[car]['orientation'] = 'H' if pos[0][0] == pos[1][0] else 'V'
        cars[car]['pos'].sort()
    return cars

def estimate_moves_to_exit(board, cars):
    # Find red car
    red_car = cars['A']
    red_row = red_car['pos'][0][0]
    red_end = red_car['pos'][-1][1]
    
    # Count blocking cars
    blocking = 0
    for j in range(red_end + 1, 6):
        if board[red_row][j] != '.' and board[red_row][j] != 'x':
            blocking += 1
    return blocking * 2 + (5 - red_end)

def get_next_states(board, cars, path_cost, bound):
    states = []
    red_car = cars['A']
    red_row = red_car['pos'][0][0]
    red_end = red_car['pos'][-1][1]
    
    def try_move(car_name, car_data, direction):
        pos = car_data['pos']
        is_horizontal = car_data['orientation'] == 'H'
        
        # Check if move is valid
        if is_horizontal:
            row = pos[0][0]
            if direction == 1 and pos[-1][1] < 5 and board[row][pos[-1][1] + 1] == '.':
                return True
            if direction == -1 and pos[0][1] > 0 and board[row][pos[0][1] - 1] == '.':
                return True
        else:
            col = pos[0][1]
            if direction == 1 and pos[-1][0] < 5 and board[pos[-1][0] + 1][col] == '.':
                return True
            if direction == -1 and pos[0][0] > 0 and board[pos[0][0] - 1][col] == '.':
                return True
        return False
    
    def make_move(board, car_name, car_data, direction):
        new_board = [list(row) for row in board]
        pos = car_data['pos']
        is_horizontal = car_data['orientation'] == 'H'
        
        # Clear current positions
        for p in pos:
            new_board[p[0]][p[1]] = '.'
        
        # Add new positions
        new_pos = []
        for p in pos:
            if is_horizontal:
                new_p = (p[0], p[1] + direction)
            else:
                new_p = (p[0] + direction, p[1])
            new_pos.append(new_p)
            new_board[new_p[0]][new_p[1]] = car_name
        
        return [''.join(row) for row in new_board]
    
    # First try moving red car if possible
    for direction in [1, -1]:
        if try_move('A', red_car, direction):
            new_board = make_move(board, 'A', red_car, direction)
            new_cost = path_cost + 1
            h_cost = estimate_moves_to_exit(new_board, cars)
            if new_cost + h_cost <= bound:
                states.append(('A', direction, new_board, new_cost))
    
    # Then try moving blocking cars
    for car_name, car_data in cars.items():
        if car_name == 'A':
            continue
        
        # Check if car is blocking
        is_blocking = any(p[0] == red_row and p[1] > red_end for p in car_data['pos'])
        
        for direction in [1, -1]:
            if try_move(car_name, car_data, direction):
                new_board = make_move(board, car_name, car_data, direction)
                new_cost = path_cost + 1
                h_cost = estimate_moves_to_exit(new_board, cars)
                if new_cost + h_cost <= bound:
                    priority = 1 if is_blocking else 2
                    states.append((car_name, direction, new_board, new_cost))
    
    return sorted(states, key=lambda x: (x[3], x[0] != 'A'))

def ida_star_search(initial_board):
    cars = get_car_data(initial_board)
    bound = estimate_moves_to_exit(initial_board, cars)
    
    def search(board, path, cost, bound):
        h_cost = estimate_moves_to_exit(board, cars)
        f_cost = cost + h_cost
        
        if f_cost > bound:
            return None, f_cost
        
        if any('A' in row and row.rindex('A') == 5 for row in board):
            return path, cost
        
        min_cost = float('inf')
        for car_name, direction, new_board, new_cost in get_next_states(board, cars, cost, bound):
            board_str = '\n'.join(new_board)
            if board_str not in visited:
                visited.add(board_str)
                move = f"{car_name}{'+' if direction > 0 else '-'}1"
                result, new_bound = search(new_board, path + [move], new_cost, bound)
                if result is not None:
                    return result, new_bound
                min_cost = min(min_cost, new_bound)
                visited.remove(board_str)
        
        return None, min_cost
    
    while True:
        visited = set()
        result, new_bound = search(initial_board, [], 0, bound)
        if result is not None:
            return result
        if new_bound == float('inf'):
            return None
        bound = new_bound

initial_board = [
    "..HBBB",
    "..HICC",
    ".AAIJK",
    ".GDDJK",
    ".GEEJL",
    "FF...L"
]

solution = ida_star_search(initial_board)
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")