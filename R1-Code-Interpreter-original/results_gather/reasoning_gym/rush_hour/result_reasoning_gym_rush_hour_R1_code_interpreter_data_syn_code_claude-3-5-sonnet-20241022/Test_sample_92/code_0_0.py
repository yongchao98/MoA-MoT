from collections import deque
import copy

def create_board():
    return [
        ['.', '.', 'x', 'B', 'B', 'B'],
        ['G', 'H', 'I', '.', 'C', 'C'],
        ['G', 'H', 'I', 'A', 'A', 'K'],
        ['D', 'D', 'D', 'J', '.', 'K'],
        ['.', '.', '.', 'J', 'E', 'E'],
        ['.', 'F', 'F', 'F', '.', 'x']
    ]

def get_vehicle_positions(board):
    vehicles = {}
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] not in '.x':
                if board[i][j] not in vehicles:
                    vehicles[board[i][j]] = []
                vehicles[board[i][j]].append((i, j))
    return vehicles

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def can_move(board, positions, di, dj):
    if is_horizontal(positions):
        if dj > 0:  # Moving right
            new_j = positions[-1][1] + dj
            return new_j < len(board[0]) and all(board[positions[0][0]][j] == '.' for j in range(positions[-1][1] + 1, new_j + 1))
        else:  # Moving left
            new_j = positions[0][1] + dj
            return new_j >= 0 and all(board[positions[0][0]][j] == '.' for j in range(new_j, positions[0][1]))
    else:  # Vertical
        if di > 0:  # Moving down
            new_i = positions[-1][0] + di
            return new_i < len(board) and all(board[i][positions[0][1]] == '.' for i in range(positions[-1][0] + 1, new_i + 1))
        else:  # Moving up
            new_i = positions[0][0] + di
            return new_i >= 0 and all(board[i][positions[0][1]] == '.' for i in range(new_i, positions[0][0]))

def move_vehicle(board, vehicle, positions, di, dj):
    new_board = [row[:] for row in board]
    # Clear old positions
    for i, j in positions:
        new_board[i][j] = '.'
    # Set new positions
    for i, j in positions:
        new_board[i + di][j + dj] = vehicle
    return new_board

def solve_puzzle():
    initial_board = create_board()
    queue = deque([(initial_board, [])])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        vehicles = get_vehicle_positions(current_board)
        
        # Check if solved (AA can move to exit)
        aa_pos = vehicles['A']
        if aa_pos[0][0] == 2 and can_move(current_board, aa_pos, 0, 4 - aa_pos[-1][1]):
            return moves + [f"A+{4 - aa_pos[-1][1]}"]
        
        # Try moving each vehicle
        for vehicle, positions in vehicles.items():
            is_hor = is_horizontal(positions)
            
            # Try moving horizontally or vertically
            for delta in [-1, 1]:
                di, dj = (0, delta) if is_hor else (delta, 0)
                if can_move(current_board, positions, di, dj):
                    new_board = move_vehicle(current_board, vehicle, positions, di, dj)
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        move = f"{vehicle}{'+' if delta > 0 else '-'}{abs(delta)}"
                        queue.append((new_board, moves + [move]))

solution = solve_puzzle()
print('<<<' + ' '.join(solution) + '>>>')