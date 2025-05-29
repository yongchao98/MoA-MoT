from collections import deque

def create_board():
    return [
        list('BBH.K.'),
        list('.GH.KL'),
        list('.GAAKL'),
        list('CCI..x'),
        list('F.IJDD'),
        list('FEEJ.x')
    ]

def get_vehicle_positions(board):
    vehicles = {}
    for i in range(6):
        for j in range(6):
            if board[i][j].isalpha():
                if board[i][j] not in vehicles:
                    vehicles[board[i][j]] = []
                vehicles[board[i][j]].append((i, j))
    return vehicles

def is_vertical(positions):
    return positions[0][1] == positions[1][1]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def can_move(board, positions, direction):
    if is_vertical(positions):  # Vertical vehicle
        col = positions[0][1]
        if direction < 0:  # Moving up
            row = min(p[0] for p in positions) - 1
            return row >= 0 and board[row][col] == '.'
        else:  # Moving down
            row = max(p[0] for p in positions) + 1
            return row < 6 and board[row][col] == '.'
    else:  # Horizontal vehicle
        row = positions[0][0]
        if direction < 0:  # Moving left
            col = min(p[1] for p in positions) - 1
            return col >= 0 and board[row][col] == '.'
        else:  # Moving right
            col = max(p[1] for p in positions) + 1
            return col < 6 and board[row][col] == '.'

def move_vehicle(board, vehicle, positions, direction):
    new_board = [row[:] for row in board]
    # Clear current positions
    for i, j in positions:
        new_board[i][j] = '.'
    
    # Set new positions
    if is_vertical(positions):
        col = positions[0][1]
        for i, j in positions:
            new_board[i + direction][j] = vehicle
    else:
        row = positions[0][0]
        for i, j in positions:
            new_board[i][j + direction] = vehicle
    
    return new_board

def solve_puzzle():
    initial_board = create_board()
    queue = deque([(initial_board, [])])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if solved
        if current_board[2][4] == 'A':
            return moves
        
        vehicles = get_vehicle_positions(current_board)
        for vehicle, positions in vehicles.items():
            positions = sorted(positions)
            
            # Try both directions
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board = move_vehicle(current_board, vehicle, positions, direction)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in visited:
                        visited.add(board_str)
                        move = f"{vehicle}{'+' if direction > 0 else '-'}{abs(direction)}"
                        queue.append((new_board, moves + [move]))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")