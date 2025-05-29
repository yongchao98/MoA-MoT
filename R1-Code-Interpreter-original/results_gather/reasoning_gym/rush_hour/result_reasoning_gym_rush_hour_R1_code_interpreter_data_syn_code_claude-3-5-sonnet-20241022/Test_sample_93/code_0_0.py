from collections import deque
import copy

def parse_board(board_str):
    return [list(row) for row in board_str.split('\n') if row]

def get_vehicle_info(board):
    vehicles = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] not in '.x':
                if board[i][j] not in vehicles:
                    vehicles[board[i][j]] = {'pos': [], 'orientation': None}
                vehicles[board[i][j]]['pos'].append((i, j))
    
    for v in vehicles:
        if vehicles[v]['pos'][0][0] == vehicles[v]['pos'][1][0]:
            vehicles[v]['orientation'] = 'H'  # horizontal
        else:
            vehicles[v]['orientation'] = 'V'  # vertical
    return vehicles

def can_move(board, vehicle, positions, direction):
    if direction in ['L', 'R'] and vehicle['orientation'] != 'H':
        return False
    if direction in ['U', 'D'] and vehicle['orientation'] != 'V':
        return False
    
    for pos in positions:
        new_pos = list(pos)
        if direction == 'L':
            new_pos[1] -= 1
        elif direction == 'R':
            new_pos[1] += 1
        elif direction == 'U':
            new_pos[0] -= 1
        else:  # D
            new_pos[0] += 1
            
        if not (0 <= new_pos[0] < len(board) and 0 <= new_pos[1] < len(board[0])):
            return False
        if board[new_pos[0]][new_pos[1]] not in ['.', pos[2]]:
            return False
    return True

def move_vehicle(board, vehicle_id, positions, direction):
    new_board = [row[:] for row in board]
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    new_positions = []
    for pos in positions:
        new_pos = list(pos)
        if direction == 'L':
            new_pos[1] -= 1
        elif direction == 'R':
            new_pos[1] += 1
        elif direction == 'U':
            new_pos[0] -= 1
        else:  # D
            new_pos[0] += 1
        new_board[new_pos[0]][new_pos[1]] = vehicle_id
        new_positions.append(tuple(new_pos))
    return new_board, new_positions

def board_to_string(board):
    return '\n'.join(''.join(row) for row in board)

def solve_puzzle(initial_board_str):
    board = parse_board(initial_board_str)
    vehicles = get_vehicle_info(board)
    
    start_state = (board_to_string(board), [])
    visited = {start_state[0]}
    queue = deque([start_state])
    
    while queue:
        current_board_str, moves = queue.popleft()
        current_board = parse_board(current_board_str)
        
        # Check if red car (AA) has reached the exit
        for i in range(len(current_board)):
            if current_board[i].count('A') == 2 and current_board[i][-1] == 'A':
                return moves
        
        # Try moving each vehicle
        for vehicle_id, vehicle in vehicles.items():
            positions = [(i, j, vehicle_id) for i in range(len(current_board)) 
                        for j in range(len(current_board[0])) 
                        if current_board[i][j] == vehicle_id]
            
            for direction, delta in [('L', -1), ('R', 1), ('U', -1), ('D', 1)]:
                if can_move(current_board, vehicle, positions, direction):
                    new_board, new_positions = move_vehicle(current_board, vehicle_id, positions, direction)
                    new_board_str = board_to_string(new_board)
                    
                    if new_board_str not in visited:
                        visited.add(new_board_str)
                        move = f"{vehicle_id}{'+' if direction in 'RD' else '-'}{1}"
                        queue.append((new_board_str, moves + [move]))
    
    return None

# Initial board
initial_board = """.BBICC
xDDI.J
x.HAAJ
G.H..K
GEE..K
...FFF"""

solution = solve_puzzle(initial_board)
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")