from collections import deque
import copy

def get_board_state():
    board = [
        ['.','.','.','.','E','.'],
        ['.','.','C','E','.','.'],
        ['A','A','C','E','.','.'],
        ['.','.','.','.','.','.'],
        ['.','.','D','.','.','.'],
        ['.','.','D','B','B','x']
    ]
    return board

def get_car_positions(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] not in ['.', 'x']:
                if board[i][j] not in cars:
                    cars[board[i][j]] = [(i,j)]
                else:
                    cars[board[i][j]].append((i,j))
    return cars

def is_vertical(positions):
    return positions[0][1] == positions[1][1]

def can_move(board, car_pos, direction):
    if len(car_pos) == 0:
        return False
    
    is_vert = is_vertical(car_pos) if len(car_pos) > 1 else None
    
    if direction in ['up', 'down'] and not is_vert:
        return False
    if direction in ['left', 'right'] and is_vert:
        return False
        
    if direction == 'up' and car_pos[0][0] > 0:
        return board[car_pos[0][0]-1][car_pos[0][1]] == '.'
    elif direction == 'down' and car_pos[-1][0] < 5:
        return board[car_pos[-1][0]+1][car_pos[-1][1]] == '.'
    elif direction == 'left' and car_pos[0][1] > 0:
        return board[car_pos[0][0]][car_pos[0][1]-1] == '.'
    elif direction == 'right' and car_pos[-1][1] < 5:
        return board[car_pos[-1][0]][car_pos[-1][1]+1] not in ['x', '#']
    return False

def move_car(board, car, direction):
    new_board = [row[:] for row in board]
    cars = get_car_positions(board)
    car_pos = cars[car]
    
    # Clear current positions
    for i, j in car_pos:
        new_board[i][j] = '.'
        
    # Calculate new positions
    new_pos = []
    for i, j in car_pos:
        if direction == 'up':
            new_pos.append((i-1, j))
        elif direction == 'down':
            new_pos.append((i+1, j))
        elif direction == 'left':
            new_pos.append((i, j-1))
        elif direction == 'right':
            new_pos.append((i, j+1))
            
    # Place car in new positions
    for i, j in new_pos:
        new_board[i][j] = car
        
    return new_board

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_puzzle():
    initial_board = get_board_state()
    queue = deque([(initial_board, [])])
    seen = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        cars = get_car_positions(current_board)
        
        # Check if solved (red car 'A' reaches right edge)
        if any(j == 5 for i, j in cars['A']):
            return moves
        
        # Try all possible moves for each car
        for car in cars:
            car_pos = cars[car]
            for direction in ['up', 'down', 'left', 'right']:
                if can_move(current_board, car_pos, direction):
                    new_board = move_car(current_board, car, direction)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in seen:
                        seen.add(board_str)
                        move_notation = f"{car}{'+' if direction in ['right', 'down'] else '-'}1"
                        queue.append((new_board, moves + [move_notation]))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")