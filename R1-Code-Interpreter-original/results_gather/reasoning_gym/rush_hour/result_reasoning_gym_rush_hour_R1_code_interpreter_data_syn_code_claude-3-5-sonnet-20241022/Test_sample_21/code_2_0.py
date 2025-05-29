from collections import deque

def get_car_positions(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j].isalpha() and board[i][j] != 'x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = []
                cars[board[i][j]].append((i, j))
    return cars

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def get_valid_moves(board, cars):
    moves = []
    # Prioritize moves that might help clear the path for car A
    priority_cars = ['A', 'E', 'G', 'H', 'B', 'D', 'C', 'F']
    
    for car in priority_cars:
        if car not in cars:
            continue
        positions = cars[car]
        is_horiz = is_horizontal(positions)
        
        if is_horiz:
            # Try moving left
            min_col = min(pos[1] for pos in positions)
            if min_col > 0 and board[positions[0][0]][min_col - 1] == '.':
                moves.append((car, -1))
            
            # Try moving right
            max_col = max(pos[1] for pos in positions)
            if max_col < 5 and board[positions[0][0]][max_col + 1] == '.':
                moves.append((car, 1))
        else:
            # Try moving up
            min_row = min(pos[0] for pos in positions)
            if min_row > 0 and board[min_row - 1][positions[0][1]] == '.':
                moves.append((car, -1))
            
            # Try moving down
            max_row = max(pos[0] for pos in positions)
            if max_row < 5 and board[max_row + 1][positions[0][1]] == '.':
                moves.append((car, 1))
    
    return moves

def make_move(board, car, direction):
    board = [list(row) for row in board]
    positions = []
    for i in range(6):
        for j in range(6):
            if board[i][j] == car:
                positions.append((i, j))
    
    is_horiz = is_horizontal(positions)
    for pos in positions:
        board[pos[0]][pos[1]] = '.'
    
    for pos in positions:
        new_row = pos[0] + (0 if is_horiz else direction)
        new_col = pos[1] + (direction if is_horiz else 0)
        board[new_row][new_col] = car
    
    return [''.join(row) for row in board]

def heuristic(board, cars):
    if 'A' not in cars:
        return float('inf')
    # Count blocking vehicles
    a_row = cars['A'][0][0]
    a_right = max(pos[1] for pos in cars['A'])
    blocking = 0
    for col in range(a_right + 1, 6):
        if board[a_row][col] != '.' and board[a_row][col] != 'x':
            blocking += 1
    return blocking

def solve_puzzle():
    initial_board = [
        '...F..',
        '..xFG.',
        'AAEFG.',
        '..E.H.',
        '.DBBH.',
        '.DCCx.'
    ]
    
    queue = [(0, initial_board, [])]  # (heuristic value, board, moves)
    visited = {''.join(initial_board)}
    
    while queue:
        queue.sort()  # Sort by heuristic value
        _, current_board, moves = queue.pop(0)
        cars = get_car_positions(current_board)
        
        # Check if solved
        if 'A' in cars and max(pos[1] for pos in cars['A']) == 5:
            return moves
        
        valid_moves = get_valid_moves(current_board, cars)
        for car, direction in valid_moves:
            new_board = make_move(current_board, car, direction)
            board_key = ''.join(new_board)
            
            if board_key not in visited:
                visited.add(board_key)
                new_moves = moves + [f"{car}{'+' if direction > 0 else ''}{direction}"]
                new_cars = get_car_positions(new_board)
                h_value = heuristic(new_board, new_cars)
                queue.append((h_value, new_board, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")