from heapq import heappush, heappop

def create_board():
    return [
        ['G', 'B', 'B', 'B', 'I', '.'],
        ['G', 'C', 'C', '.', 'I', '.'],
        ['A', 'A', 'H', '.', 'I', 'J'],
        ['.', 'x', 'H', 'D', 'D', 'J'],
        ['E', 'E', 'E', '.', '.', 'J'],
        ['.', '.', 'F', 'F', 'F', '.']
    ]

def get_car_positions(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] != '.' and board[i][j] != 'x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = [(i, j)]
                else:
                    cars[board[i][j]].append((i, j))
    return cars

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def better_heuristic(board, cars):
    # Focus on red car's path
    red_row = cars['A'][0][0]
    red_right = max(pos[1] for pos in cars['A'])
    
    # Count direct blockers
    direct_blockers = 0
    blocking_cars = set()
    for col in range(red_right + 1, 6):
        if board[red_row][col] != '.' and board[red_row][col] != 'x':
            direct_blockers += 3
            blocking_cars.add(board[red_row][col])
    
    # Add penalty for cars blocking the blockers
    secondary_blockers = 0
    for car in blocking_cars:
        positions = cars[car]
        if not is_horizontal(positions):  # vertical car
            col = positions[0][1]
            for pos in positions:
                row = pos[0]
                # Check if this blocker is blocked
                if row > 0 and board[row-1][col] != '.':
                    secondary_blockers += 1
                if row < 5 and board[row+1][col] != '.':
                    secondary_blockers += 1
    
    return direct_blockers + secondary_blockers

def can_move(board, car_pos, direction):
    if is_horizontal(car_pos):
        row = car_pos[0][0]
        if direction > 0:
            col = max(pos[1] for pos in car_pos) + 1
            return col < 6 and board[row][col] == '.'
        else:
            col = min(pos[1] for pos in car_pos) - 1
            return col >= 0 and board[row][col] == '.'
    else:
        col = car_pos[0][1]
        if direction > 0:
            row = max(pos[0] for pos in car_pos) + 1
            return row < 6 and board[row][col] == '.'
        else:
            row = min(pos[0] for pos in car_pos) - 1
            return row >= 0 and board[row][col] == '.'

def move_car(board, car, car_pos, direction):
    new_board = [row[:] for row in board]
    for pos in car_pos:
        new_board[pos[0]][pos[1]] = '.'
    if is_horizontal(car_pos):
        row = car_pos[0][0]
        for pos in car_pos:
            new_board[row][pos[1] + direction] = car
    else:
        col = car_pos[0][1]
        for pos in car_pos:
            new_board[pos[0] + direction][col] = car
    return new_board

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_puzzle():
    initial_board = create_board()
    initial_cars = get_car_positions(initial_board)
    pq = [(better_heuristic(initial_board, initial_cars), 0, initial_board, [])]
    seen = {board_to_string(initial_board)}
    
    while pq:
        _, cost, current_board, moves = heappop(pq)
        cars = get_car_positions(current_board)
        
        # Check if solved
        if any(pos[1] == 4 for pos in cars['A']):
            return moves
        
        # Prioritize moves that directly affect red car's path
        red_row = cars['A'][0][0]
        red_right = max(pos[1] for pos in cars['A'])
        blocking_cols = set(j for j in range(red_right + 1, 6) 
                          if current_board[red_row][j] != '.')
        
        for car, positions in cars.items():
            # Prioritize moves for cars blocking the red car
            is_blocker = any(pos[1] in blocking_cols and pos[0] == red_row 
                           for pos in positions)
            priority = 0 if is_blocker else 1
            
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board = move_car(current_board, car, positions, direction)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in seen:
                        seen.add(board_str)
                        new_cars = get_car_positions(new_board)
                        h_score = better_heuristic(new_board, new_cars)
                        move_str = f"{car}{'+' if direction > 0 else '-'}1"
                        heappush(pq, (cost + 1 + h_score + priority, 
                                    cost + 1, new_board, moves + [move_str]))
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")