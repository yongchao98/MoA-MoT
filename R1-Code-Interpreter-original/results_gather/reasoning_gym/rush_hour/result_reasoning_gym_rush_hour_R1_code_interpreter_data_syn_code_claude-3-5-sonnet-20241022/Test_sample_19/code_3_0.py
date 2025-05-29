from collections import deque

def get_vehicle_positions(board):
    vehicles = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] != '.' and board[i][j] != 'x':
                if board[i][j] not in vehicles:
                    vehicles[board[i][j]] = [(i, j)]
                else:
                    vehicles[board[i][j]].append((i, j))
    return vehicles

def is_horizontal(positions):
    return len(positions) == 1 or positions[0][0] == positions[1][0]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def can_move(board, positions, direction):
    if is_horizontal(positions):
        if direction > 0:  # right
            last_col = max(p[1] for p in positions)
            if last_col + 1 >= len(board[0]):
                return False
            return board[positions[0][0]][last_col + 1] == '.'
        else:  # left
            first_col = min(p[1] for p in positions)
            if first_col - 1 < 0:
                return False
            return board[positions[0][0]][first_col - 1] == '.'
    else:  # vertical
        if direction > 0:  # down
            last_row = max(p[0] for p in positions)
            if last_row + 1 >= len(board):
                return False
            return board[last_row + 1][positions[0][1]] == '.'
        else:  # up
            first_row = min(p[0] for p in positions)
            if first_row - 1 < 0:
                return False
            return board[first_row - 1][positions[0][1]] == '.'

def move_vehicle(board, vehicle, positions, direction):
    new_board = [row[:] for row in board]
    # Remove old positions
    for row, col in positions:
        new_board[row][col] = '.'
    
    # Add new positions
    if is_horizontal(positions):
        for row, col in positions:
            new_board[row][col + direction] = vehicle
    else:
        for row, col in positions:
            new_board[row + direction][col] = vehicle
    
    return new_board

def solve_puzzle():
    initial_board = [
        ['.', 'F', 'B', 'B', 'J', 'K'],
        ['.', 'F', '.', 'I', 'J', 'K'],
        ['x', 'A', 'A', 'I', '.', 'L'],
        ['C', 'C', 'H', 'I', '.', 'L'],
        ['E', 'G', 'H', '.', '.', '.'],
        ['E', 'G', 'D', 'D', 'D', '.']
    ]
    
    queue = deque([(initial_board, [])])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        vehicles = get_vehicle_positions(current_board)
        
        # Check if solved (AA reaches the right edge)
        if 'A' in vehicles:
            aa_positions = vehicles['A']
            if max(pos[1] for pos in aa_positions) == len(current_board[0]) - 1:
                return moves
        
        # Try all possible moves for each vehicle
        for vehicle, positions in vehicles.items():
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board = move_vehicle(current_board, vehicle, positions, direction)
                    board_string = board_to_string(new_board)
                    
                    if board_string not in visited:
                        visited.add(board_string)
                        new_moves = moves + [f"{vehicle}{'+' if direction > 0 else '-'}1"]
                        queue.append((new_board, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")