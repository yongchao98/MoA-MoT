from heapq import heappush, heappop

def get_car_info(board):
    cars = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            car = board[i][j]
            if car != '.' and car != 'x':
                if car not in cars:
                    cars[car] = {'positions': set([(i, j)]), 'length': 1}
                else:
                    cars[car]['positions'].add((i, j))
                    cars[car]['length'] += 1
    return cars

def get_heuristic(board, car_info):
    # Calculate minimum moves needed to get red car to exit
    red_positions = car_info['A']['positions']
    red_row = list(red_positions)[0][0]
    red_rightmost = max(x[1] for x in red_positions)
    blocking_count = 0
    
    # Count blocking cars and their positions
    for j in range(red_rightmost + 1, len(board[0])):
        if board[red_row][j] != '.' and board[red_row][j] != 'A':
            blocking_count += 1
    
    return blocking_count

def get_moves(board, car_info):
    moves = []
    for car, info in car_info.items():
        positions = list(info['positions'])
        is_horizontal = positions[0][0] == positions[1][0] if len(positions) > 1 else None
        
        if is_horizontal:
            # Try horizontal moves
            row = positions[0][0]
            leftmost = min(p[1] for p in positions)
            rightmost = max(p[1] for p in positions)
            
            # Try move left
            if leftmost > 0 and board[row][leftmost-1] == '.':
                moves.append((car, -1))
            # Try move right
            if rightmost < len(board[0])-1 and board[row][rightmost+1] == '.':
                moves.append((car, 1))
        else:
            # Try vertical moves
            col = positions[0][1]
            topmost = min(p[0] for p in positions)
            bottommost = max(p[0] for p in positions)
            
            # Try move up
            if topmost > 0 and board[topmost-1][col] == '.':
                moves.append((car, -1))
            # Try move down
            if bottommost < len(board)-1 and board[bottommost+1][col] == '.':
                moves.append((car, 1))
    
    return moves

def apply_move(board, car_info, move):
    car, direction = move
    positions = list(car_info[car]['positions'])
    new_board = [list(row) for row in board]
    
    # Clear current positions
    for i, j in positions:
        new_board[i][j] = '.'
    
    # Set new positions
    is_horizontal = positions[0][0] == positions[1][0] if len(positions) > 1 else None
    new_positions = []
    for i, j in positions:
        new_i = i + (0 if is_horizontal else direction)
        new_j = j + (direction if is_horizontal else 0)
        new_board[new_i][new_j] = car
        new_positions.append((new_i, new_j))
    
    return [''.join(row) for row in new_board]

def solve_puzzle():
    initial_board = [
        "BBBCC.",
        "..JKDD",
        "AAJKLM",
        "I.EELM",
        "IFF..N",
        "GGHHHN"
    ]
    
    visited = set()
    pq = [(0, initial_board, [])]  # (priority, board, moves)
    visited.add(tuple(initial_board))
    
    while pq:
        _, current_board, moves = heappop(pq)
        car_info = get_car_info(current_board)
        
        # Check if solved
        red_positions = car_info['A']['positions']
        if max(x[1] for x in red_positions) == len(current_board[0]) - 1:
            return moves
        
        # Get possible moves
        possible_moves = get_moves(current_board, car_info)
        
        # Try each move
        for move in possible_moves:
            new_board = apply_move(current_board, car_info, move)
            board_tuple = tuple(new_board)
            
            if board_tuple not in visited:
                visited.add(board_tuple)
                car, direction = move
                move_str = f"{car}{'+' if direction > 0 else ''}{direction}"
                new_moves = moves + [move_str]
                
                # Calculate priority based on heuristic
                priority = len(new_moves) + get_heuristic(new_board, get_car_info(new_board))
                heappush(pq, (priority, new_board, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")