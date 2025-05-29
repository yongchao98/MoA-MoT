from collections import deque
import copy

def get_car_positions(board):
    cars = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] not in ['.', 'x']:
                if board[i][j] not in cars:
                    cars[board[i][j]] = [(i, j)]
                else:
                    cars[board[i][j]].append((i, j))
    return cars

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def can_move(board, car, car_positions, delta_row, delta_col):
    for pos in car_positions:
        new_row, new_col = pos[0] + delta_row, pos[1] + delta_col
        if not (0 <= new_row < 6 and 0 <= new_col < 6):
            return False
        if board[new_row][new_col] not in ['.', car] or board[new_row][new_col] == 'x':
            return False
    return True

def move_car(board, car, positions, delta_row, delta_col):
    new_board = [list(row) for row in board]
    # Clear old positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    # Set new positions
    for pos in positions:
        new_board[pos[0] + delta_row][pos[1] + delta_col] = car
    return new_board

def solve_puzzle():
    initial_board = [
        ['B', 'B', 'I', 'C', 'C', 'x'],
        ['.', '.', 'I', 'J', 'D', 'D'],
        ['A', 'A', 'I', 'J', '.', 'K'],
        ['H', 'E', 'E', 'F', 'F', 'K'],
        ['H', '.', 'G', 'G', '.', '.'],
        ['.', '.', '.', '.', '.', '.']
    ]
    
    queue = deque([(initial_board, [])])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if red car (AA) can reach exit
        cars = get_car_positions(current_board)
        if 'A' in cars:
            red_car = cars['A']
            # Check if path to exit is clear
            if red_car[0][0] == 2 and all(current_board[2][red_car[1][1]+1:5] == ['.' for _ in range(red_car[1][1]+1, 5)]):
                return moves + [f"A+{5-red_car[1][1]}"]
        
        for car, positions in cars.items():
            is_hor = is_horizontal(positions)
            
            if is_hor:
                # Try moving left
                if can_move(current_board, car, positions, 0, -1):
                    new_board = move_car(current_board, car, positions, 0, -1)
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        queue.append((new_board, moves + [f"{car}-1"]))
                
                # Try moving right
                if can_move(current_board, car, positions, 0, 1):
                    new_board = move_car(current_board, car, positions, 0, 1)
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        queue.append((new_board, moves + [f"{car}+1"]))
            else:
                # Try moving up
                if can_move(current_board, car, positions, -1, 0):
                    new_board = move_car(current_board, car, positions, -1, 0)
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        queue.append((new_board, moves + [f"{car}-1"]))
                
                # Try moving down
                if can_move(current_board, car, positions, 1, 0):
                    new_board = move_car(current_board, car, positions, 1, 0)
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        queue.append((new_board, moves + [f"{car}+1"]))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")