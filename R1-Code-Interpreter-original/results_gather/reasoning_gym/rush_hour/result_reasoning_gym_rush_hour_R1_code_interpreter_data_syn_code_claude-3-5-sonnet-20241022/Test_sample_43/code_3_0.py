def get_car_positions(board):
    cars = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] != '.' and board[i][j] != 'x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = [(i, j)]
                else:
                    cars[board[i][j]].append((i, j))
    return cars

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def try_move(board, car, positions, direction, horizontal):
    new_board = [list(row) for row in board]
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Try new positions
    try:
        for pos in positions:
            new_pos = (pos[0], pos[1] + direction) if horizontal else (pos[0] + direction, pos[1])
            if new_pos[0] < 0 or new_pos[1] < 0:
                return None
            if new_board[new_pos[0]][new_pos[1]] not in ['.']:
                return None
            new_board[new_pos[0]][new_pos[1]] = car
    except IndexError:
        return None
    
    return [''.join(row) for row in new_board]

def solve_puzzle(board, moves=None, seen=None):
    if moves is None:
        moves = []
    if seen is None:
        seen = {tuple(board)}
    
    cars = get_car_positions(board)
    
    # Check if solved
    aa_positions = cars['A']
    if max(p[1] for p in aa_positions) == len(board[0])-2:
        return moves
    
    # Try moving cars that are blocking AA first
    aa_row = aa_positions[0][0]
    blocking_cars = set()
    for j in range(max(p[1] for p in aa_positions) + 1, len(board[0])):
        if board[aa_row][j] not in '.x':
            blocking_cars.add(board[aa_row][j])
    
    # First try moving blocking cars
    for car in blocking_cars:
        positions = cars[car]
        horizontal = is_horizontal(positions)
        
        for direction in [-1, 1]:
            new_board = try_move(board, car, positions, direction, horizontal)
            if new_board and tuple(new_board) not in seen:
                seen.add(tuple(new_board))
                move = f"{car}{'+' if direction > 0 else ''}{direction}"
                result = solve_puzzle(new_board, moves + [move], seen)
                if result:
                    return result
    
    # Then try moving AA
    positions = cars['A']
    for direction in [-1, 1]:
        new_board = try_move(board, 'A', positions, direction, True)
        if new_board and tuple(new_board) not in seen:
            seen.add(tuple(new_board))
            move = f"A{'+' if direction > 0 else ''}{direction}"
            result = solve_puzzle(new_board, moves + [move], seen)
            if result:
                return result
    
    # Finally try other cars
    for car, positions in cars.items():
        if car not in blocking_cars and car != 'A':
            horizontal = is_horizontal(positions)
            for direction in [-1, 1]:
                new_board = try_move(board, car, positions, direction, horizontal)
                if new_board and tuple(new_board) not in seen:
                    seen.add(tuple(new_board))
                    move = f"{car}{'+' if direction > 0 else ''}{direction}"
                    result = solve_puzzle(new_board, moves + [move], seen)
                    if result:
                        return result
    
    return None

initial_board = [
    "BBB.xK",
    "CCDDDK",
    "HAAJ..",
    "H..JEE",
    "..IFFx",
    "GGI..."
]

solution = solve_puzzle(initial_board)
if solution:
    print('<<<' + ' '.join(solution) + '>>>')
else:
    print("No solution found")