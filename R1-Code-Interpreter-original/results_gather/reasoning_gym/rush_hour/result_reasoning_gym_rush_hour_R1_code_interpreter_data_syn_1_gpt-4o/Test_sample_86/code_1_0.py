from collections import deque

def is_valid_move(board, car, direction, steps):
    car_positions = [(r, c) for r, row in enumerate(board) for c, cell in enumerate(row) if cell == car]
    if not car_positions:
        return False
    
    if direction in ['right', 'left']:
        row = car_positions[0][0]
        if direction == 'right':
            new_positions = [(row, c + steps) for _, c in car_positions]
        else:
            new_positions = [(row, c - steps) for _, c in car_positions]
    else:
        col = car_positions[0][1]
        if direction == 'down':
            new_positions = [(r + steps, col) for r, _ in car_positions]
        else:
            new_positions = [(r - steps, col) for r, _ in car_positions]
    
    for r, c in new_positions:
        if r < 0 or r >= len(board) or c < 0 or c >= len(board[0]) or (board[r][c] != '.' and board[r][c] != car):
            return False
    return True

def move_car(board, car, direction, steps):
    new_board = [list(row) for row in board]
    car_positions = [(r, c) for r, row in enumerate(board) for c, cell in enumerate(row) if cell == car]
    
    if direction in ['right', 'left']:
        row = car_positions[0][0]
        if direction == 'right':
            for _, c in sorted(car_positions, reverse=True):
                new_board[row][c + steps] = car
                new_board[row][c] = '.'
        else:
            for _, c in sorted(car_positions):
                new_board[row][c - steps] = car
                new_board[row][c] = '.'
    else:
        col = car_positions[0][1]
        if direction == 'down':
            for r, _ in sorted(car_positions, reverse=True):
                new_board[r + steps][col] = car
                new_board[r][col] = '.'
        else:
            for r, _ in sorted(car_positions):
                new_board[r - steps][col] = car
                new_board[r][col] = '.'
    
    return [''.join(row) for row in new_board]

def solve_puzzle(initial_board):
    queue = deque([(initial_board, [])])
    visited = set()
    visited.add(tuple(initial_board))
    
    while queue:
        board, moves = queue.popleft()
        
        if board[2][5] == 'A':
            return moves
        
        cars = set(cell for row in board for cell in row if cell not in '.x')
        for car in cars:
            for direction in ['right', 'left', 'up', 'down']:
                for steps in range(1, 3):  # Limit steps to 2 for efficiency
                    if is_valid_move(board, car, direction, steps):
                        new_board = move_car(board, car, direction, steps)
                        new_state = tuple(new_board)
                        if new_state not in visited:
                            visited.add(new_state)
                            new_moves = moves + [f"{car}{'+' if direction in ['right', 'down'] else '-'}{steps}"]
                            queue.append((new_board, new_moves))
    
    return None

initial_board = [
    ".BB.xM",
    "..JCCM",
    "AAJ.LM",
    "DD.KL.",
    "IEEKFF",
    "IGGGHH"
]

solution = solve_puzzle(initial_board)
print(solution)