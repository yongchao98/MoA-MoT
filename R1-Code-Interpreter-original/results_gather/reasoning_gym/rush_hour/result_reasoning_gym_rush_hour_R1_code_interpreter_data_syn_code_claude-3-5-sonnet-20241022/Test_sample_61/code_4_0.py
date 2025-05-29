from collections import deque

# Initial board state
initial_board = [
    list("BBBCCM"),
    list("DDEEEM"),
    list("AAJK.."),
    list("I.JKFF"),
    list("I.GGLx"),
    list("HH..L.")
]

def get_car_positions(board):
    cars = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] not in ['.', 'x']:
                if board[i][j] not in cars:
                    cars[board[i][j]] = []
                cars[board[i][j]].append((i, j))
    return cars

def get_valid_moves(board):
    moves = []
    cars = get_car_positions(board)
    
    for car, positions in cars.items():
        # Determine if car is horizontal or vertical
        is_horizontal = positions[0][0] == positions[-1][0]
        
        if is_horizontal:
            row = positions[0][0]
            left = min(p[1] for p in positions)
            right = max(p[1] for p in positions)
            
            # Check left move
            if left > 0 and board[row][left-1] == '.':
                moves.append((car, -1))
            # Check right move
            if right < len(board[0])-1 and board[row][right+1] == '.':
                moves.append((car, 1))
        else:
            col = positions[0][1]
            top = min(p[0] for p in positions)
            bottom = max(p[0] for p in positions)
            
            # Check up move
            if top > 0 and board[top-1][col] == '.':
                moves.append((car, -1))
            # Check down move
            if bottom < len(board)-1 and board[bottom+1][col] == '.':
                moves.append((car, 1))
    
    return moves

def apply_move(board, move):
    car, direction = move
    new_board = [row[:] for row in board]
    cars = get_car_positions(board)
    positions = cars[car]
    is_horizontal = positions[0][0] == positions[-1][0]
    
    # Remove car from old position
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Add car to new position
    for pos in positions:
        new_row = pos[0]
        new_col = pos[1]
        if is_horizontal:
            new_col += direction
        else:
            new_row += direction
        new_board[new_row][new_col] = car
    
    return new_board

def is_solved(board):
    cars = get_car_positions(board)
    if 'A' in cars:
        red_car_positions = cars['A']
        return max(pos[1] for pos in red_car_positions) == 4
    return False

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve():
    start = (initial_board, [])
    queue = deque([start])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        if is_solved(current_board):
            return moves
        
        for move in get_valid_moves(current_board):
            new_board = apply_move(current_board, move)
            board_str = board_to_string(new_board)
            
            if board_str not in visited:
                visited.add(board_str)
                new_moves = moves + [move]
                queue.append((new_board, new_moves))
    
    return None

# Execute solution
solution = solve()
if solution:
    formatted_solution = ''
    for car, direction in solution:
        sign = '+' if direction > 0 else '-'
        formatted_solution += f"{car}{sign}1 "
    print(f"<<<{formatted_solution.strip()}>>>")
else:
    print("<<<No solution found>>>")