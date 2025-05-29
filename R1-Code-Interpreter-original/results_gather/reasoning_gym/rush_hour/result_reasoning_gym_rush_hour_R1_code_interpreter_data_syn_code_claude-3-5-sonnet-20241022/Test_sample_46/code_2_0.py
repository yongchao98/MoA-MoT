from collections import deque

def get_car_positions(board):
    cars = {}
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j].isalpha() and board[i][j] != 'x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = []
                cars[board[i][j]].append((i, j))
    return {car: sorted(positions) for car, positions in cars.items()}

def is_horizontal(positions):
    return positions[0][0] == positions[-1][0]

def get_blocking_cars(board):
    cars = get_car_positions(board)
    blocking = set()
    
    # Find red car (AA) position
    red_car = cars['A']
    red_row = red_car[0][0]
    red_right = red_car[-1][1]
    
    # Check all cars blocking the path to exit
    for car, positions in cars.items():
        if car != 'A':
            if not is_horizontal(positions):  # vertical car
                col = positions[0][1]
                if col > red_right and any(pos[0] == red_row for pos in positions):
                    blocking.add(car)
    
    return blocking

def solve_puzzle(board):
    def apply_move(board, car, direction):
        new_board = [list(row) for row in board]
        positions = get_car_positions(board)[car]
        
        if is_horizontal(positions):
            row = positions[0][0]
            if direction < 0:  # left
                new_board[row][positions[0][1]-1] = car
                new_board[row][positions[-1][1]] = '.'
            else:  # right
                new_board[row][positions[-1][1]+1] = car
                new_board[row][positions[0][1]] = '.'
        else:  # vertical
            col = positions[0][1]
            if direction < 0:  # up
                new_board[positions[0][0]-1][col] = car
                new_board[positions[-1][0]][col] = '.'
            else:  # down
                new_board[positions[-1][0]+1][col] = car
                new_board[positions[0][0]][col] = '.'
        
        return [''.join(row) for row in new_board]

    def get_valid_moves(board):
        moves = []
        cars = get_car_positions(board)
        blocking = get_blocking_cars(board)
        
        # Prioritize moving blocking cars
        for car in cars:
            positions = cars[car]
            if car in blocking or car == 'A':
                if is_horizontal(positions):
                    row = positions[0][0]
                    # Try left
                    if positions[0][1] > 0 and board[row][positions[0][1]-1] == '.':
                        moves.append((car, -1))
                    # Try right
                    if positions[-1][1] < len(board[0])-1 and board[row][positions[-1][1]+1] == '.':
                        moves.append((car, 1))
                else:
                    col = positions[0][1]
                    # Try up
                    if positions[0][0] > 0 and board[positions[0][0]-1][col] == '.':
                        moves.append((car, -1))
                    # Try down
                    if positions[-1][0] < len(board)-1 and board[positions[-1][0]+1][col] == '.':
                        moves.append((car, 1))
        
        return moves

    def is_solved(board):
        for row in board:
            if 'A' in row and row.index('A') == len(row)-2:
                return True
        return False

    queue = deque([(board, [])])
    seen = {tuple(board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        if is_solved(current_board):
            return moves
        
        for car, direction in get_valid_moves(current_board):
            new_board = apply_move(current_board, car, direction)
            board_tuple = tuple(new_board)
            
            if board_tuple not in seen:
                seen.add(board_tuple)
                new_moves = moves + [f"{car}{'+' if direction > 0 else ''}{direction}"]
                queue.append((new_board, new_moves))
    
    return None

# Initial board
board = [
    "xBBB.I",
    ".CC..I",
    "..GAAI",
    "..GHDD",
    "..GHEE",
    ".FFx.."
]

solution = solve_puzzle(board)
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")