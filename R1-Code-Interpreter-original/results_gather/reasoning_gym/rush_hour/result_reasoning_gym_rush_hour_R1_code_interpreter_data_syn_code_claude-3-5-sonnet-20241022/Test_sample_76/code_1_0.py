from collections import deque

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

def get_car_orientation(positions):
    # Returns 'H' for horizontal, 'V' for vertical
    return 'H' if positions[0][0] == positions[1][0] else 'V'

def can_move(board, positions, direction):
    orientation = get_car_orientation(positions)
    if orientation == 'H':
        row = positions[0][0]
        if direction > 0:  # right
            max_col = max(p[1] for p in positions)
            return max_col + 1 < len(board[0]) and board[row][max_col + 1] == '.'
        else:  # left
            min_col = min(p[1] for p in positions)
            return min_col - 1 >= 0 and board[row][min_col - 1] == '.'
    else:  # vertical
        col = positions[0][1]
        if direction > 0:  # down
            max_row = max(p[0] for p in positions)
            return max_row + 1 < len(board) and board[max_row + 1][col] == '.'
        else:  # up
            min_row = min(p[0] for p in positions)
            return min_row - 1 >= 0 and board[min_row - 1][col] == '.'

def move_car(board, car, positions, direction):
    new_board = [list(row) for row in board]
    orientation = get_car_orientation(positions)
    
    # Clear current positions
    for i, j in positions:
        new_board[i][j] = '.'
    
    # Add new positions
    if orientation == 'H':
        row = positions[0][0]
        min_col = min(p[1] for p in positions)
        for i in range(len(positions)):
            new_board[row][min_col + i + direction] = car
    else:
        col = positions[0][1]
        min_row = min(p[0] for p in positions)
        for i in range(len(positions)):
            new_board[min_row + i + direction][col] = car
    
    return [''.join(row) for row in new_board]

def get_possible_moves(board, cars):
    moves = []
    for car, positions in cars.items():
        for direction in [-1, 1]:
            if can_move(board, positions, direction):
                new_board = move_car(board, car, positions, direction)
                moves.append((car, direction, new_board))
    return moves

def solve_puzzle():
    initial_board = [
        "FBBBCC",
        "F.G.DD",
        "AAG.HI",
        "..G.HI",
        ".....J",
        "EEE..J"
    ]
    
    queue = deque([(initial_board, [])])
    visited = {''.join(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        cars = get_car_positions(current_board)
        
        # Check if solved
        if 'A' in cars:
            red_car_pos = cars['A']
            if max(p[1] for p in red_car_pos) == len(current_board[0]) - 1:
                return moves
        
        # Get all possible moves
        possible_moves = get_possible_moves(current_board, cars)
        
        # Try each move
        for car, direction, new_board in possible_moves:
            board_str = ''.join(new_board)
            if board_str not in visited:
                visited.add(board_str)
                new_moves = moves + [f"{car}{'+' if direction > 0 else '-'}1"]
                queue.append((new_board, new_moves))
                
                # If this is a move of car 'A', prioritize it by adding it to the front
                if car == 'A':
                    queue.appendleft((new_board, new_moves))
    
    return None

# Solve and print the solution
solution = solve_puzzle()
if solution:
    print('<<<' + ' '.join(solution) + '>>>')
else:
    print("No solution found")