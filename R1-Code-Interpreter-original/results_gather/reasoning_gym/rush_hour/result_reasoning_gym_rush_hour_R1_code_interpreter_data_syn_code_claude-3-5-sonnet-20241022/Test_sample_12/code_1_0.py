def get_car_info(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] != '.' and board[i][j] != 'x':
                car = board[i][j]
                if car not in cars:
                    cars[car] = {'positions': [(i, j)], 'length': 1}
                else:
                    cars[car]['positions'].append((i, j))
                    cars[car]['length'] += 1
    
    for car in cars:
        pos = cars[car]['positions']
        pos.sort()
        cars[car]['orientation'] = 'H' if pos[0][0] == pos[-1][0] else 'V'
        cars[car]['row'] = pos[0][0]
        cars[car]['col'] = pos[0][1]
    return cars

def is_valid_move(board, car_info, distance):
    positions = car_info['positions']
    orientation = car_info['orientation']
    length = car_info['length']
    
    if orientation == 'H':
        row = positions[0][0]
        col = positions[0][1]
        if distance > 0:  # right
            return col + length - 1 + distance < 6 and all(board[row][col + length + i] == '.' for i in range(distance))
        else:  # left
            return col + distance >= 0 and all(board[row][col + i] == '.' for i in range(distance, 0))
    else:  # vertical
        row = positions[0][0]
        col = positions[0][1]
        if distance > 0:  # down
            return row + length - 1 + distance < 6 and all(board[row + length + i][col] == '.' for i in range(distance))
        else:  # up
            return row + distance >= 0 and all(board[row + i][col] == '.' for i in range(distance, 0))

def make_move(board, car, car_info, distance):
    new_board = [list(row) for row in board]
    positions = car_info['positions']
    orientation = car_info['orientation']
    length = car_info['length']
    
    # Clear old positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Set new positions
    if orientation == 'H':
        row = positions[0][0]
        col = positions[0][1]
        for i in range(length):
            new_board[row][col + distance + i] = car
    else:
        row = positions[0][0]
        col = positions[0][1]
        for i in range(length):
            new_board[row + distance + i][col] = car
    
    return [''.join(row) for row in new_board]

def solve_puzzle(board):
    cars = get_car_info(board)
    moves = []
    
    # First, clear the path for car A
    # Move H up
    if 'H' in cars:
        moves.append('H-2')
        board = make_move(board, 'H', cars['H'], -2)
        cars = get_car_info(board)
    
    # Move C left
    if 'C' in cars:
        moves.append('C-2')
        board = make_move(board, 'C', cars['C'], -2)
        cars = get_car_info(board)
    
    # Move I left
    if 'I' in cars:
        moves.append('I-2')
        board = make_move(board, 'I', cars['I'], -2)
        cars = get_car_info(board)
    
    # Move B left
    if 'B' in cars:
        moves.append('B-2')
        board = make_move(board, 'B', cars['B'], -2)
        cars = get_car_info(board)
    
    # Move G left
    if 'G' in cars:
        moves.append('G-2')
        board = make_move(board, 'G', cars['G'], -2)
        cars = get_car_info(board)
    
    # Finally move A to exit
    moves.append('A+4')
    
    return moves

# Initial board
board = [
    "GBBIJK",
    "G..IJK",
    "AAHI..",
    "..HCCC",
    "..xDD.",
    "EEEFF."
]

solution = solve_puzzle(board)
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")