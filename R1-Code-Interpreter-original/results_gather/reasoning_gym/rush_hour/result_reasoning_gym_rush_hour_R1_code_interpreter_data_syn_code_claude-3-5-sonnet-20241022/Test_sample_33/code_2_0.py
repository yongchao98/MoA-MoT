from collections import deque

def get_car_info(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] not in '.x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = {'positions': [(i, j)], 'orientation': None}
                else:
                    cars[board[i][j]]['positions'].append((i, j))
                    # Determine orientation
                    if cars[board[i][j]]['positions'][0][0] == i:
                        cars[board[i][j]]['orientation'] = 'horizontal'
                    else:
                        cars[board[i][j]]['orientation'] = 'vertical'
    return cars

def can_move(board, car_info, car, direction):
    positions = car_info[car]['positions']
    orientation = car_info[car]['orientation']
    
    if orientation == 'horizontal':
        row = positions[0][0]
        if direction == -1:  # Left
            col = min(p[1] for p in positions)
            return col > 0 and board[row][col-1] == '.'
        else:  # Right
            col = max(p[1] for p in positions)
            return col < 5 and board[row][col+1] == '.'
    else:  # vertical
        col = positions[0][1]
        if direction == -1:  # Up
            row = min(p[0] for p in positions)
            return row > 0 and board[row-1][col] == '.'
        else:  # Down
            row = max(p[0] for p in positions)
            return row < 5 and board[row+1][col] == '.'

def make_move(board, car_info, car, direction):
    new_board = [list(row) for row in board]
    positions = car_info[car]['positions']
    orientation = car_info[car]['orientation']
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Calculate new positions
    new_positions = []
    for pos in positions:
        if orientation == 'horizontal':
            new_pos = (pos[0], pos[1] + direction)
        else:
            new_pos = (pos[0] + direction, pos[1])
        new_board[new_pos[0]][new_pos[1]] = car
        new_positions.append(new_pos)
    
    new_car_info = dict(car_info)
    new_car_info[car]['positions'] = new_positions
    
    return [''.join(row) for row in new_board], new_car_info

def solve_rush_hour(initial_board):
    board = initial_board.split('\n')
    car_info = get_car_info(board)
    
    # BFS with state pruning
    visited = {''.join(board)}
    queue = deque([(board, car_info, [])])
    
    while queue:
        current_board, current_car_info, moves = queue.popleft()
        
        # Check if solved (AA reaches exit)
        aa_positions = current_car_info['A']['positions']
        if max(p[1] for p in aa_positions) == 5:
            return moves
        
        # Try moving each car
        for car in current_car_info:
            for direction in [-1, 1]:
                if can_move(current_board, current_car_info, car, direction):
                    new_board, new_car_info = make_move(current_board, current_car_info, car, direction)
                    board_str = ''.join(new_board)
                    
                    if board_str not in visited:
                        visited.add(board_str)
                        new_moves = moves + [(car, direction)]
                        queue.append((new_board, new_car_info, new_moves))
    
    return None

# Initial board
initial_board = """..xBBK
GCCJ.K
GAAJ..
DDIJ..
.HIEE.
.HFFx."""

# Solve and format output
solution = solve_rush_hour(initial_board)
if solution:
    move_str = ''
    for car, spaces in solution:
        move_str += f'{car}{"+"+str(spaces) if spaces > 0 else str(spaces)} '
    print(f"<<<{move_str.strip()}>>>")
else:
    print("No solution found")