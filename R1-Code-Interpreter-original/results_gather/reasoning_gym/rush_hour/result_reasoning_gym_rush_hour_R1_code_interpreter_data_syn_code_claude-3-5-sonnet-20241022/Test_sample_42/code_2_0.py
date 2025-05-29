from heapq import heappush, heappop

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

def get_valid_moves(board, cars):
    moves = []
    for car, positions in cars.items():
        if len(positions) < 2:
            continue
            
        horizontal = is_horizontal(positions)
        if horizontal:
            # Try moving left
            leftmost = min(p[1] for p in positions)
            if leftmost > 0 and board[positions[0][0]][leftmost-1] == '.':
                moves.append((car, -1))
            # Try moving right
            rightmost = max(p[1] for p in positions)
            if rightmost < 5 and board[positions[0][0]][rightmost+1] == '.':
                moves.append((car, 1))
        else:
            # Try moving up
            topmost = min(p[0] for p in positions)
            if topmost > 0 and board[topmost-1][positions[0][1]] == '.':
                moves.append((car, -1))
            # Try moving down
            bottommost = max(p[0] for p in positions)
            if bottommost < 5 and board[bottommost+1][positions[0][1]] == '.':
                moves.append((car, 1))
    return moves

def apply_move(board, cars, move):
    car, direction = move
    positions = cars[car]
    new_board = [list(row) for row in board]
    horizontal = is_horizontal(positions)
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Add new positions
    new_positions = []
    for pos in positions:
        if horizontal:
            new_pos = (pos[0], pos[1] + direction)
        else:
            new_pos = (pos[0] + direction, pos[1])
        new_board[new_pos[0]][new_pos[1]] = car
        new_positions.append(new_pos)
    
    return [''.join(row) for row in new_board]

def manhattan_distance(cars, target='A'):
    if target in cars:
        rightmost = max(p[1] for p in cars[target])
        return 4 - rightmost
    return float('inf')

def solve_puzzle(initial_board):
    pq = [(0, 0, initial_board, [])]
    seen = {tuple(initial_board)}
    
    while pq:
        _, moves_count, current_board, moves = heappop(pq)
        cars = get_car_positions(current_board)
        
        # Check if solved
        if 'A' in cars and max(p[1] for p in cars['A']) == 4:
            return moves
        
        valid_moves = get_valid_moves(current_board, cars)
        for move in valid_moves:
            new_board = apply_move(current_board, cars, move)
            board_tuple = tuple(new_board)
            
            if board_tuple not in seen:
                seen.add(board_tuple)
                car, direction = move
                move_str = f"{car}{'+' if direction > 0 else ''}{direction}"
                new_cars = get_car_positions(new_board)
                priority = manhattan_distance(new_cars) + moves_count + 1
                heappush(pq, (priority, moves_count + 1, new_board, moves + [move_str]))
    
    return None

# Initial board
initial_board = [
    "G..x..",
    "GBBJ.L",
    "AAIJ.L",
    "CCIDDL",
    ".HEEK.",
    ".HFFK."
]

solution = solve_puzzle(initial_board)
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")