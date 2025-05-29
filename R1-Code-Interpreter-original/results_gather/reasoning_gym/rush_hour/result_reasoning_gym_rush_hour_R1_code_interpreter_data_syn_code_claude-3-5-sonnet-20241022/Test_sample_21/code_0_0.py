from collections import deque
import copy

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
    for car, positions in cars.items():
        is_horiz = is_horizontal(positions)
        
        # Try moving left/up
        can_move = True
        min_pos = min(positions[0][1] if is_horiz else positions[0][0])
        if min_pos > 0:
            check_pos = (positions[0][0], min_pos - 1) if is_horiz else (min_pos - 1, positions[0][1])
            if board[check_pos[0]][check_pos[1]] == '.':
                moves.append((car, -1))
                
        # Try moving right/down
        max_pos = max(positions[-1][1] if is_horiz else positions[-1][0])
        if max_pos < 5:
            check_pos = (positions[0][0], max_pos + 1) if is_horiz else (max_pos + 1, positions[0][1])
            if board[check_pos[0]][check_pos[1]] == '.':
                moves.append((car, 1))
    
    return moves

def apply_move(board, cars, move):
    car, direction = move
    positions = cars[car]
    new_board = [list(row) for row in board]
    is_horiz = is_horizontal(positions)
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Add new positions
    new_positions = []
    for pos in positions:
        new_pos = (pos[0], pos[1] + direction) if is_horiz else (pos[0] + direction, pos[1])
        new_board[new_pos[0]][new_pos[1]] = car
        new_positions.append(new_pos)
    
    return [''.join(row) for row in new_board], new_positions

def solve_puzzle():
    initial_board = [
        '...F..',
        '..xFG.',
        'AAEFG.',
        '..E.H.',
        '.DBBH.',
        '.DCCx.'
    ]
    
    queue = deque([(initial_board, [], [])])
    visited = set(''.join(initial_board))
    
    while queue:
        current_board, moves, cars_state = queue.popleft()
        current_cars = get_car_positions(current_board)
        
        # Check if solved (AA reaches right edge)
        if 'A' in current_cars:
            a_positions = current_cars['A']
            if max(pos[1] for pos in a_positions) == 5:
                return moves
        
        valid_moves = get_valid_moves(current_board, current_cars)
        for move in valid_moves:
            new_board, new_positions = apply_move(current_board, current_cars, move)
            board_key = ''.join(new_board)
            
            if board_key not in visited:
                visited.add(board_key)
                new_moves = moves + [f"{move[0]}{'+' if move[1] > 0 else ''}{move[1]}"]
                queue.append((new_board, new_moves, new_positions))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")