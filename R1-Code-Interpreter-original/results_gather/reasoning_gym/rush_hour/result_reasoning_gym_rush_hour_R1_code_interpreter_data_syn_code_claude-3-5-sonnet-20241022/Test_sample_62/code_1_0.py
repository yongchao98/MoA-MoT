from collections import deque
import copy

def create_board():
    board = [
        ['G', 'B', 'B', '.', 'L', 'M'],
        ['G', 'H', 'I', '.', 'L', 'M'],
        ['G', 'H', 'I', 'A', 'A', 'N'],
        ['C', 'C', 'J', 'K', '.', 'N'],
        ['.', '.', 'J', 'K', 'D', 'D'],
        ['.', 'E', 'E', 'F', 'F', '.']
    ]
    return board

def find_vehicles(board):
    vehicles = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] != '.' and board[i][j] != 'x':
                if board[i][j] not in vehicles:
                    vehicles[board[i][j]] = {'pos': [(i, j)], 'orientation': None}
                else:
                    vehicles[board[i][j]]['pos'].append((i, j))
    
    for v in vehicles:
        if vehicles[v]['pos'][0][0] == vehicles[v]['pos'][1][0]:
            vehicles[v]['orientation'] = 'H'  # horizontal
        else:
            vehicles[v]['orientation'] = 'V'  # vertical
    return vehicles

def get_possible_moves(board, vehicles):
    moves = []
    for v in vehicles:
        if vehicles[v]['orientation'] == 'H':
            # Try moving left
            leftmost = min(p[1] for p in vehicles[v]['pos'])
            if leftmost > 0 and board[vehicles[v]['pos'][0][0]][leftmost-1] == '.':
                moves.append((v, -1))
            # Try moving right
            rightmost = max(p[1] for p in vehicles[v]['pos'])
            if rightmost < 5 and board[vehicles[v]['pos'][0][0]][rightmost+1] == '.':
                moves.append((v, 1))
        else:
            # Try moving up
            topmost = min(p[0] for p in vehicles[v]['pos'])
            if topmost > 0 and board[topmost-1][vehicles[v]['pos'][0][1]] == '.':
                moves.append((v, -1))
            # Try moving down
            bottommost = max(p[0] for p in vehicles[v]['pos'])
            if bottommost < 5 and board[bottommost+1][vehicles[v]['pos'][0][1]] == '.':
                moves.append((v, 1))
    return moves

def apply_move(board, vehicles, move):
    vehicle, direction = move
    new_board = [row[:] for row in board]
    new_vehicles = copy.deepcopy(vehicles)
    
    # Clear current position
    for pos in vehicles[vehicle]['pos']:
        new_board[pos[0]][pos[1]] = '.'
    
    # Update position
    if vehicles[vehicle]['orientation'] == 'H':
        new_vehicles[vehicle]['pos'] = [(p[0], p[1] + direction) for p in vehicles[vehicle]['pos']]
    else:
        new_vehicles[vehicle]['pos'] = [(p[0] + direction, p[1]) for p in vehicles[vehicle]['pos']]
    
    # Place vehicle in new position
    for pos in new_vehicles[vehicle]['pos']:
        new_board[pos[0]][pos[1]] = vehicle
    
    return new_board, new_vehicles

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_puzzle():
    initial_board = create_board()
    initial_vehicles = find_vehicles(initial_board)
    
    queue = deque([(initial_board, initial_vehicles, [])])
    seen = {board_to_string(initial_board)}
    
    while queue:
        current_board, current_vehicles, moves = queue.popleft()
        
        # Check if solved (AA is at the exit)
        if any(pos == (2, 4) for pos in current_vehicles['A']['pos']):
            return moves
        
        # Generate and try all possible moves
        possible_moves = get_possible_moves(current_board, current_vehicles)
        for move in possible_moves:
            new_board, new_vehicles = apply_move(current_board, current_vehicles, move)
            board_str = board_to_string(new_board)
            
            if board_str not in seen:
                seen.add(board_str)
                new_moves = moves + [(move[0], move[1])]
                queue.append((new_board, new_vehicles, new_moves))
    
    return None

# Solve and format the solution
solution = solve_puzzle()
if solution:
    formatted_solution = ''
    for vehicle, direction in solution:
        formatted_solution += f'{vehicle}{"+"+str(direction) if direction > 0 else str(direction)} '
    print(formatted_solution.strip())