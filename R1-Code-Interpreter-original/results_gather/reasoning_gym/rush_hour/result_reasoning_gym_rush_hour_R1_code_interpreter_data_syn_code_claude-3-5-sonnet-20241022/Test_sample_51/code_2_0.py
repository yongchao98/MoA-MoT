from collections import deque, defaultdict

def get_board_string(board):
    return ''.join(''.join(row) for row in board)

def get_car_data(board):
    cars = defaultdict(list)
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] not in ['.', 'x']:
                cars[board[i][j]].append((i, j))
    return {k: sorted(v) for k, v in cars.items()}

def is_blocking_red_car(car_pos, red_car_row):
    return any(pos[0] == red_car_row and pos[1] > max(p[1] for p in red_car_pos['A']) for pos in car_pos.values())

def get_moves(board, car_pos):
    moves = []
    for car, positions in car_pos.items():
        is_horizontal = positions[0][0] == positions[1][0]
        if is_horizontal:
            row = positions[0][0]
            left = min(p[1] for p in positions) - 1
            right = max(p[1] for p in positions) + 1
            
            if left >= 0 and board[row][left] == '.':
                moves.append((car, -1, 1 if car == 'A' else 2))
            if right < len(board[0]) and board[row][right] == '.':
                moves.append((car, 1, 0 if car == 'A' else 2))
        else:
            col = positions[0][1]
            top = min(p[0] for p in positions) - 1
            bottom = max(p[0] for p in positions) + 1
            
            if top >= 0 and board[top][col] == '.':
                moves.append((car, -1, 2))
            if bottom < len(board) and board[bottom][col] == '.':
                moves.append((car, 1, 2))
    return sorted(moves, key=lambda x: x[2])

def apply_move(board, car_pos, car, direction):
    new_board = [list(row) for row in board]
    positions = car_pos[car]
    is_horizontal = positions[0][0] == positions[1][0]
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Add new positions
    new_positions = []
    for pos in positions:
        if is_horizontal:
            new_pos = (pos[0], pos[1] + direction)
        else:
            new_pos = (pos[0] + direction, pos[1])
        new_board[new_pos[0]][new_pos[1]] = car
        new_positions.append(new_pos)
    
    return new_board

def solve_rush_hour():
    initial_board = [
        list("BBBKCC"),
        list("DDJK.L"),
        list("I.JAAL"),
        list("IEE.xM"),
        list("FF...M"),
        list("GG.xHH")
    ]
    
    visited = set()
    queue = [(initial_board, [])]
    visited.add(get_board_string(initial_board))
    
    while queue:
        current_board, moves = queue.pop(0)
        car_pos = get_car_data(current_board)
        
        # Check if solved
        if 'A' in car_pos and max(p[1] for p in car_pos['A']) == len(current_board[0]) - 1:
            return moves
        
        # Get and sort moves by priority
        valid_moves = get_moves(current_board, car_pos)
        
        # Try each move
        for car, direction, _ in valid_moves:
            new_board = apply_move(current_board, car_pos, car, direction)
            board_string = get_board_string(new_board)
            
            if board_string not in visited:
                visited.add(board_string)
                new_moves = moves + [(car, direction)]
                
                # Prioritize moves that directly help solve the puzzle
                if car == 'A' and direction > 0:
                    queue.insert(0, (new_board, new_moves))
                else:
                    queue.append((new_board, new_moves))
    
    return None

# Solve and format the output
solution = solve_rush_hour()
if solution:
    formatted_moves = []
    for car, direction in solution:
        sign = '+' if direction > 0 else '-'
        formatted_moves.append(f"{car}{sign}1")
    print("<<<" + " ".join(formatted_moves) + ">>>")
else:
    print("No solution found")