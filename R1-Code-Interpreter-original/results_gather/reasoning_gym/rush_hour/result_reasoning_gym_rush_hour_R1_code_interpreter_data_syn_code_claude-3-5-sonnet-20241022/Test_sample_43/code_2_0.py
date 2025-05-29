from heapq import heappush, heappop

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

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

def heuristic(board, cars):
    # Distance of red car from exit
    aa_positions = cars['A']
    exit_pos = (2, 5)  # Exit position
    return manhattan_distance((aa_positions[0][0], max(p[1] for p in aa_positions)), exit_pos)

def get_blocking_cars(board, cars):
    aa_positions = cars['A']
    aa_row = aa_positions[0][0]
    aa_rightmost = max(p[1] for p in aa_positions)
    blocking = set()
    
    # Check cars directly blocking AA
    for j in range(aa_rightmost + 1, len(board[0])):
        if board[aa_row][j] not in '.x':
            blocking.add(board[aa_row][j])
    return blocking

def solve_puzzle(initial_board):
    cars = get_car_positions(initial_board)
    initial_state = (tuple(initial_board), ())
    seen = {initial_state[0]}
    pq = [(heuristic(initial_board, cars), 0, initial_state)]
    
    while pq:
        _, cost, (current_board, moves) = heappop(pq)
        current_board = list(current_board)
        cars = get_car_positions(current_board)
        
        # Check if solved
        aa_positions = cars['A']
        if max(p[1] for p in aa_positions) == len(current_board[0])-2:
            return list(moves)
            
        blocking_cars = get_blocking_cars(current_board, cars)
        
        # Generate moves
        for car in cars:
            positions = cars[car]
            horizontal = is_horizontal(positions)
            
            # Prioritize blocking cars and AA
            if car not in blocking_cars and car != 'A' and len(blocking_cars) > 0:
                continue
                
            if horizontal:
                # Try moving left
                leftmost = min(p[1] for p in positions)
                if leftmost > 0 and current_board[positions[0][0]][leftmost-1] == '.':
                    new_board = [list(row) for row in current_board]
                    for pos in positions:
                        new_board[pos[0]][pos[1]] = '.'
                        new_board[pos[0]][pos[1]-1] = car
                    new_board = tuple(''.join(row) for row in new_board)
                    
                    if new_board not in seen:
                        seen.add(new_board)
                        new_moves = moves + (f"{car}-1",)
                        new_cars = get_car_positions(new_board)
                        heappush(pq, (cost + 1 + heuristic(new_board, new_cars), 
                                    cost + 1, 
                                    (new_board, new_moves)))
                
                # Try moving right
                rightmost = max(p[1] for p in positions)
                if rightmost < len(current_board[0])-1 and current_board[positions[0][0]][rightmost+1] == '.':
                    new_board = [list(row) for row in current_board]
                    for pos in positions:
                        new_board[pos[0]][pos[1]] = '.'
                        new_board[pos[0]][pos[1]+1] = car
                    new_board = tuple(''.join(row) for row in new_board)
                    
                    if new_board not in seen:
                        seen.add(new_board)
                        new_moves = moves + (f"{car}+1",)
                        new_cars = get_car_positions(new_board)
                        heappush(pq, (cost + 1 + heuristic(new_board, new_cars), 
                                    cost + 1, 
                                    (new_board, new_moves)))
            else:
                # Try moving up
                topmost = min(p[0] for p in positions)
                if topmost > 0 and current_board[topmost-1][positions[0][1]] == '.':
                    new_board = [list(row) for row in current_board]
                    for pos in positions:
                        new_board[pos[0]][pos[1]] = '.'
                        new_board[pos[0]-1][pos[1]] = car
                    new_board = tuple(''.join(row) for row in new_board)
                    
                    if new_board not in seen:
                        seen.add(new_board)
                        new_moves = moves + (f"{car}-1",)
                        new_cars = get_car_positions(new_board)
                        heappush(pq, (cost + 1 + heuristic(new_board, new_cars), 
                                    cost + 1, 
                                    (new_board, new_moves)))
                
                # Try moving down
                bottommost = max(p[0] for p in positions)
                if bottommost < len(current_board)-1 and current_board[bottommost+1][positions[0][1]] == '.':
                    new_board = [list(row) for row in current_board]
                    for pos in positions:
                        new_board[pos[0]][pos[1]] = '.'
                        new_board[pos[0]+1][pos[1]] = car
                    new_board = tuple(''.join(row) for row in new_board)
                    
                    if new_board not in seen:
                        seen.add(new_board)
                        new_moves = moves + (f"{car}+1",)
                        new_cars = get_car_positions(new_board)
                        heappush(pq, (cost + 1 + heuristic(new_board, new_cars), 
                                    cost + 1, 
                                    (new_board, new_moves)))
    
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