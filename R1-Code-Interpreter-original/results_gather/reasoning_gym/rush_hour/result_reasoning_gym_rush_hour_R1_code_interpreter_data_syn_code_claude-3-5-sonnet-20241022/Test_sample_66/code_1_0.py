def get_car_info(board):
    cars = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] != '.' and board[i][j] != 'x':
                car = board[i][j]
                if car not in cars:
                    cars[car] = {'pos': [(i, j)]}
                else:
                    cars[car]['pos'].append((i, j))
                    # Set orientation
                    cars[car]['horizontal'] = cars[car]['pos'][0][0] == cars[car]['pos'][1][0]
    return cars

def can_move(board, car_info, direction):
    positions = car_info['pos']
    if car_info['horizontal']:
        row = positions[0][0]
        if direction > 0:  # right
            col = max(pos[1] for pos in positions) + 1
            return col < len(board[0]) and board[row][col] == '.'
        else:  # left
            col = min(pos[1] for pos in positions) - 1
            return col >= 0 and board[row][col] == '.'
    else:  # vertical
        col = positions[0][1]
        if direction > 0:  # down
            row = max(pos[0] for pos in positions) + 1
            return row < len(board) and board[row][col] == '.'
        else:  # up
            row = min(pos[0] for pos in positions) - 1
            return row >= 0 and board[row][col] == '.'

def move_car(board, car, car_info, direction):
    new_board = [row[:] for row in board]
    positions = car_info['pos']
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Set new positions
    if car_info['horizontal']:
        row = positions[0][0]
        for pos in positions:
            new_board[row][pos[1] + direction] = car
    else:
        col = positions[0][1]
        for pos in positions:
            new_board[pos[0] + direction][col] = car
    
    return new_board

def is_solved(board):
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] == 'A' and j == len(board[i])-1:
                return True
    return False

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_puzzle_dfs(board, max_depth, visited, moves):
    if len(moves) > max_depth:
        return None
    
    if is_solved(board):
        return moves
    
    cars = get_car_info(board)
    
    # First priority: Try moving the red car (A) to the right if possible
    if 'A' in cars:
        if can_move(board, cars['A'], 1):  # Try moving right
            new_board = move_car(board, 'A', cars['A'], 1)
            board_str = board_to_string(new_board)
            if board_str not in visited:
                visited.add(board_str)
                result = solve_puzzle_dfs(new_board, max_depth, visited, moves + ['A+1'])
                if result:
                    return result
                visited.remove(board_str)
    
    # Second priority: Move cars blocking the red car's path
    aa_row = cars['A']['pos'][0][0]
    aa_right = max(pos[1] for pos in cars['A']['pos'])
    blocking_cars = set()
    
    for col in range(aa_right + 1, len(board[0])):
        if board[aa_row][col] not in ['.', 'x']:
            blocking_cars.add(board[aa_row][col])
    
    # Try moving blocking cars first
    for car in blocking_cars:
        for direction in [-1, 1]:
            if can_move(board, cars[car], direction):
                new_board = move_car(board, car, cars[car], direction)
                board_str = board_to_string(new_board)
                if board_str not in visited:
                    visited.add(board_str)
                    result = solve_puzzle_dfs(new_board, max_depth, visited, 
                                           moves + [f"{car}{'+' if direction > 0 else '-'}{abs(direction)}"])
                    if result:
                        return result
                    visited.remove(board_str)
    
    # Finally, try moving other cars that might help indirectly
    for car, info in cars.items():
        if car != 'A' and car not in blocking_cars:
            for direction in [-1, 1]:
                if can_move(board, info, direction):
                    new_board = move_car(board, car, info, direction)
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        result = solve_puzzle_dfs(new_board, max_depth, visited, 
                                               moves + [f"{car}{'+' if direction > 0 else '-'}{abs(direction)}"])
                        if result:
                            return result
                        visited.remove(board_str)
    
    return None

def solve_puzzle():
    initial_board = [
        ['x','.','G','B','B','J'],
        ['C','C','G','H','I','J'],
        ['F','A','A','H','I','K'],
        ['F','D','D','.','I','K'],
        ['E','E','.','x','.','.'],
        ['.','.','.','.','.','.']
    ]
    
    # Iterative deepening
    for depth in range(1, 15):
        visited = {board_to_string(initial_board)}
        solution = solve_puzzle_dfs(initial_board, depth, visited, [])
        if solution:
            return solution
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")