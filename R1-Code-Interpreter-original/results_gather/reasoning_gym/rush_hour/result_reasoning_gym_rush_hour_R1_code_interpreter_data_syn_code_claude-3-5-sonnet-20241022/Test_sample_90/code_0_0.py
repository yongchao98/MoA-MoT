from collections import deque
import copy

def create_board():
    board = [
        ['B', 'B', 'B', 'K', 'C', 'C'],
        ['D', 'D', 'J', 'K', '.', 'L'],
        ['I', '.', 'J', 'A', 'A', 'L'],
        ['I', 'E', 'E', '.', 'x', 'M'],
        ['F', 'F', '.', '.', '.', 'M'],
        ['G', 'G', 'H', 'H', '.', 'x']
    ]
    return board

def get_vehicle_info(board):
    vehicles = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] not in ['.', 'x']:
                if board[i][j] not in vehicles:
                    vehicles[board[i][j]] = {'positions': [(i, j)], 'orientation': None}
                else:
                    vehicles[board[i][j]]['positions'].append((i, j))
    
    for v in vehicles:
        pos = vehicles[v]['positions']
        vehicles[v]['orientation'] = 'H' if pos[0][0] == pos[-1][0] else 'V'
    return vehicles

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def can_move(board, vehicle, positions, direction):
    new_positions = []
    for i, j in positions:
        ni, nj = i + direction[0], j + direction[1]
        if not (0 <= ni < 6 and 0 <= nj < 6) or board[ni][nj] not in ['.', vehicle]:
            return False, []
        new_positions.append((ni, nj))
    return True, new_positions

def solve_puzzle():
    initial_board = create_board()
    vehicles = get_vehicle_info(initial_board)
    
    queue = deque([(initial_board, [])])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if solved (AA is at the exit)
        if current_board[2][4] == 'A' and current_board[2][5] == 'A':
            return moves
        
        # Try moving each vehicle
        for vehicle, info in vehicles.items():
            directions = [(0, 1), (0, -1)] if info['orientation'] == 'H' else [(1, 0), (-1, 0)]
            
            for direction in directions:
                can_move_result, new_positions = can_move(current_board, vehicle, info['positions'], direction)
                if can_move_result:
                    new_board = [row[:] for row in current_board]
                    
                    # Clear old positions
                    for i, j in info['positions']:
                        new_board[i][j] = '.'
                    
                    # Set new positions
                    for i, j in new_positions:
                        new_board[i][j] = vehicle
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        move_str = f"{vehicle}{'+' if direction[0] + direction[1] > 0 else '-'}1"
                        queue.append((new_board, moves + [move_str]))

    return None

solution = solve_puzzle()
if solution:
    print('<<<' + ' '.join(solution) + '>>>')
else:
    print("No solution found")