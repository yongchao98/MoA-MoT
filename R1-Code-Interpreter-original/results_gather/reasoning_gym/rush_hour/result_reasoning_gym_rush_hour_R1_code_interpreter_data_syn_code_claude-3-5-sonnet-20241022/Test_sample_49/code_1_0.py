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

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def get_valid_moves(board, cars):
    moves = []
    # Prioritize cars I and J as they directly block car A
    priority_cars = ['I', 'J']
    remaining_cars = [car for car in cars.keys() if car not in priority_cars]
    
    for car in priority_cars + remaining_cars:
        positions = cars[car]
        horizontal = is_horizontal(positions)
        
        if horizontal:
            row = positions[0][0]
            left = min(p[1] for p in positions) - 1
            right = max(p[1] for p in positions) + 1
            
            if left >= 0 and board[row][left] == '.':
                moves.append((car, -1))
            if right < len(board[0]) and board[row][right] == '.':
                moves.append((car, 1))
        else:
            col = positions[0][1]
            top = min(p[0] for p in positions) - 1
            bottom = max(p[0] for p in positions) + 1
            
            if top >= 0 and board[top][col] == '.':
                moves.append((car, -1))
            if bottom < len(board) and board[bottom][col] == '.':
                moves.append((car, 1))
    
    return moves

def apply_move(board, cars, car, direction):
    new_board = [list(row) for row in board]
    positions = cars[car]
    horizontal = is_horizontal(positions)
    
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    new_positions = []
    for pos in positions:
        if horizontal:
            new_pos = (pos[0], pos[1] + direction)
        else:
            new_pos = (pos[0] + direction, pos[1])
        new_board[new_pos[0]][new_pos[1]] = car
        new_positions.append(new_pos)
    
    return [''.join(row) for row in new_board], new_positions

def solve_puzzle(initial_board):
    board = initial_board.split('\n')
    cars = get_car_positions(board)
    visited = set()
    queue = deque([(board, cars, [])])
    
    while queue:
        current_board, current_cars, moves = queue.popleft()
        board_str = '\n'.join(current_board)
        
        if board_str in visited:
            continue
        visited.add(board_str)
        
        if any(pos[1] == len(board[0])-1 for pos in current_cars['A']):
            return moves
        
        for car, direction in get_valid_moves(current_board, current_cars):
            new_board, new_positions = apply_move(current_board, current_cars, car, direction)
            new_cars = dict(current_cars)
            new_cars[car] = new_positions
            
            move_str = f"{car}{'+' if direction > 0 else ''}{direction}"
            queue.append((new_board, new_cars, moves + [move_str]))
    
    return None

initial_board = """.....x
BB.xJ.
FAAIJ.
F.HICC
G.HDDK
GEEE.K"""

solution = solve_puzzle(initial_board)
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")