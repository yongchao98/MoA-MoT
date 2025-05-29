from collections import deque
import copy

def print_board(board):
    for row in board:
        print(''.join(row))

def move_vehicle(board, vehicle, direction):
    # Find vehicle positions
    positions = []
    for i in range(6):
        for j in range(6):
            if board[i][j] == vehicle:
                positions.append((i, j))
    
    # Determine orientation (horizontal or vertical)
    is_horizontal = positions[0][0] == positions[1][0]
    
    # Check if move is valid
    if is_horizontal:
        if direction > 0:  # Move right
            if positions[-1][1] + direction >= 6:
                return None
            for pos in positions:
                for d in range(1, direction + 1):
                    if pos[1] + d >= 6 or board[pos[0]][pos[1] + d] not in ['.']:
                        return None
        else:  # Move left
            if positions[0][1] + direction < 0:
                return None
            for pos in positions:
                for d in range(direction, 0):
                    if pos[1] + d < 0 or board[pos[0]][pos[1] + d] not in ['.']:
                        return None
    else:  # Vertical
        if direction > 0:  # Move down
            if positions[-1][0] + direction >= 6:
                return None
            for pos in positions:
                for d in range(1, direction + 1):
                    if pos[0] + d >= 6 or board[pos[0] + d][pos[1]] not in ['.']:
                        return None
        else:  # Move up
            if positions[0][0] + direction < 0:
                return None
            for pos in positions:
                for d in range(direction, 0):
                    if pos[0] + d < 0 or board[pos[0] + d][pos[1]] not in ['.']:
                        return None
    
    # Make move
    new_board = [row[:] for row in board]
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    for pos in positions:
        if is_horizontal:
            new_board[pos[0]][pos[1] + direction] = vehicle
        else:
            new_board[pos[0] + direction][pos[1]] = vehicle
    
    return new_board

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
        
        # Check if solved (AA is at the exit)
        if any(current_board[2][i:i+2] == ['A', 'A'] for i in [4, 5]):
            return moves
        
        # Try all possible moves for all vehicles
        for vehicle in set(''.join(''.join(row) for row in current_board)):
            if vehicle in '.x':
                continue
                
            # Try moving in both directions
            for direction in [-1, 1]:
                # Try different distances
                for distance in [1, 2, 3]:
                    new_board = move_vehicle(current_board, vehicle, direction * distance)
                    if new_board:
                        board_str = board_to_string(new_board)
                        if board_str not in seen:
                            seen.add(board_str)
                            new_moves = moves + [f"{vehicle}{'+' if direction > 0 else '-'}{distance}"]
                            queue.append((new_board, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")