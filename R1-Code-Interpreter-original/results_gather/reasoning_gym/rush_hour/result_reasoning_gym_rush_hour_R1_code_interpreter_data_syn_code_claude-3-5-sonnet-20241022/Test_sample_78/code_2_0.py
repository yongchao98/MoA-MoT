from collections import deque

def print_board(board):
    for row in board:
        print(''.join(row))

def get_vehicle_positions(board):
    vehicles = {}
    for i in range(6):
        for j in range(6):
            if board[i][j] not in '.x':
                if board[i][j] not in vehicles:
                    vehicles[board[i][j]] = []
                vehicles[board[i][j]].append((i, j))
    return vehicles

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def make_move(board, vehicle, positions, delta):
    new_board = [row[:] for row in board]
    horizontal = is_horizontal(positions)
    
    # Clear current positions
    for i, j in positions:
        new_board[i][j] = '.'
    
    # Set new positions
    for i, j in positions:
        if horizontal:
            new_board[i][j + delta] = vehicle
        else:
            new_board[i + delta][j] = vehicle
            
    return new_board

def can_move(board, positions, delta):
    horizontal = is_horizontal(positions)
    if horizontal:
        row = positions[0][0]
        if delta > 0:  # Moving right
            rightmost = max(p[1] for p in positions)
            if rightmost + delta >= 6:
                return False
            return all(board[row][p[1] + delta] == '.' for p in positions)
        else:  # Moving left
            leftmost = min(p[1] for p in positions)
            if leftmost + delta < 0:
                return False
            return all(board[row][p[1] + delta] == '.' for p in positions)
    else:  # Vertical
        col = positions[0][1]
        if delta > 0:  # Moving down
            bottommost = max(p[0] for p in positions)
            if bottommost + delta >= 6:
                return False
            return all(board[p[0] + delta][col] == '.' for p in positions)
        else:  # Moving up
            topmost = min(p[0] for p in positions)
            if topmost + delta < 0:
                return False
            return all(board[p[0] + delta][col] == '.' for p in positions)

def is_solved(board):
    # Check if red car (AA) is at the exit
    return board[2][4:6] == ['A', 'A']

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_puzzle():
    initial_board = [
        list("BBBH.."),
        list("F..H.."),
        list("FAAI.J"),
        list("CCGI.J"),
        list("..GDDJ"),
        list("xxEE..")
    ]
    
    queue = deque([(initial_board, [])])
    seen = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        if is_solved(current_board):
            return moves
        
        vehicles = get_vehicle_positions(current_board)
        
        for vehicle, positions in vehicles.items():
            # Sort positions for consistency
            positions.sort()
            
            # Determine possible move ranges
            max_delta = 3 if is_horizontal(positions) else 3
            
            for delta in range(-max_delta, max_delta + 1):
                if delta == 0:
                    continue
                
                if can_move(current_board, positions, delta):
                    new_board = make_move(current_board, vehicle, positions, delta)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in seen:
                        seen.add(board_str)
                        new_moves = moves + [f"{vehicle}{'+' if delta > 0 else '-'}{abs(delta)}"]
                        queue.append((new_board, new_moves))
                        
                        # Debug: Print when finding new valid moves
                        if vehicle == 'A':
                            print(f"Found move: {new_moves[-1]}")
                            print_board(new_board)
                            print()
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")