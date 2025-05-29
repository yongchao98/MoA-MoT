from collections import deque
import copy

def print_board(board):
    for row in board:
        print(''.join(row))

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
    for car, positions in cars.items():
        is_horiz = is_horizontal(positions)
        
        if is_horiz:
            # Try moving left
            leftmost = min(p[1] for p in positions)
            if leftmost > 0 and board[positions[0][0]][leftmost-1] == '.':
                moves.append((car, -1))
            
            # Try moving right
            rightmost = max(p[1] for p in positions)
            if rightmost < len(board[0])-1 and board[positions[0][0]][rightmost+1] == '.':
                moves.append((car, 1))
        else:
            # Try moving up
            topmost = min(p[0] for p in positions)
            if topmost > 0 and board[topmost-1][positions[0][1]] == '.':
                moves.append((car, -1))
            
            # Try moving down
            bottommost = max(p[0] for p in positions)
            if bottommost < len(board)-1 and board[bottommost+1][positions[0][1]] == '.':
                moves.append((car, 1))
    
    return moves

def apply_move(board, cars, move):
    car, direction = move
    positions = cars[car]
    new_board = [list(row) for row in board]
    new_cars = copy.deepcopy(cars)
    
    is_horiz = is_horizontal(positions)
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    new_positions = []
    for pos in positions:
        if is_horiz:
            new_pos = (pos[0], pos[1] + direction)
        else:
            new_pos = (pos[0] + direction, pos[1])
        new_positions.append(new_pos)
        new_board[new_pos[0]][new_pos[1]] = car
    
    new_cars[car] = new_positions
    return [''.join(row) for row in new_board], new_cars

def solve_puzzle():
    initial_board = [
        "BBBCC.",
        "..JKDD",
        "AAJKLM",
        "I.EELM",
        "IFF..N",
        "GGHHHN"
    ]
    
    visited = set()
    queue = deque([(initial_board, [], {})])
    visited.add(tuple(initial_board))
    
    while queue:
        current_board, moves, _ = queue.popleft()
        cars = get_car_positions(current_board)
        
        # Check if solved (red car 'A' reaches exit)
        if any(pos for pos in cars['A'] if pos[1] == len(current_board[0])-1):
            return moves
        
        for move in get_valid_moves(current_board, cars):
            new_board, new_cars = apply_move(current_board, cars, move)
            board_tuple = tuple(new_board)
            
            if board_tuple not in visited:
                visited.add(board_tuple)
                car, direction = move
                move_str = f"{car}{'+' if direction > 0 else ''}{direction}"
                queue.append((new_board, moves + [move_str], new_cars))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")