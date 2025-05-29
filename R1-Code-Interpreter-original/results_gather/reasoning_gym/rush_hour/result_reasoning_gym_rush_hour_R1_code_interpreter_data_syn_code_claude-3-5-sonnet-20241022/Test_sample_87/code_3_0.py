from collections import deque
import heapq

def parse_board(board_str):
    return [list(row) for row in board_str.split('\n') if row]

def get_vehicle_positions(board):
    vehicles = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] != '.' and board[i][j] != 'x':
                if board[i][j] not in vehicles:
                    vehicles[board[i][j]] = []
                vehicles[board[i][j]].append((i, j))
    return {k: sorted(v) for k, v in vehicles.items()}

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def get_blocking_vehicles(board, red_car_row, red_car_right):
    blocking = []
    for col in range(red_car_right + 1, len(board[0])):
        if board[red_car_row][col] != '.' and board[red_car_row][col] != 'A':
            blocking.append(board[red_car_row][col])
    return blocking

def can_move(board, positions, direction):
    is_horiz = is_horizontal(positions)
    if is_horiz:
        row = positions[0][0]
        if direction > 0:
            col = max(pos[1] for pos in positions) + 1
            return col < len(board[0]) and board[row][col] == '.'
        else:
            col = min(pos[1] for pos in positions) - 1
            return col >= 0 and board[row][col] == '.'
    else:
        col = positions[0][1]
        if direction > 0:
            row = max(pos[0] for pos in positions) + 1
            return row < len(board) and board[row][col] == '.'
        else:
            row = min(pos[0] for pos in positions) - 1
            return row >= 0 and board[row][col] == '.'

def move_vehicle(board, vehicle, positions, direction):
    new_board = [row[:] for row in board]
    is_horiz = is_horizontal(positions)
    
    for row, col in positions:
        new_board[row][col] = '.'
    
    for row, col in positions:
        new_row = row + (0 if is_horiz else direction)
        new_col = col + (direction if is_horiz else 0)
        new_board[new_row][new_col] = vehicle
    
    return new_board

def board_to_string(board):
    return '\n'.join(''.join(row) for row in board)

def solve_puzzle():
    initial_board = parse_board("""G..IBB
G.HI.J
G.HAAJ
CCDD.K
EEEFFK
......""")
    
    # Priority queue elements: (priority, board, moves)
    queue = [(0, initial_board, [])]
    seen = {board_to_string(initial_board)}
    
    while queue:
        _, current_board, moves = heapq.heappop(queue)
        vehicles = get_vehicle_positions(current_board)
        
        # Get red car position
        red_car = vehicles['A']
        red_row = red_car[0][0]
        red_right = max(pos[1] for pos in red_car)
        
        # Check if solved
        if red_right == len(current_board[0]) - 1:
            return moves
        
        # Get blocking vehicles
        blocking = get_blocking_vehicles(current_board, red_row, red_right)
        
        # Try moves for each vehicle
        for vehicle, positions in vehicles.items():
            for direction in [-1, 1]:
                if can_move(current_board, positions, direction):
                    new_board = move_vehicle(current_board, vehicle, positions, direction)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in seen:
                        seen.add(board_str)
                        new_moves = moves + [f"{vehicle}{'+' if direction > 0 else '-'}{abs(direction)}"]
                        
                        # Calculate priority
                        priority = len(new_moves)
                        if vehicle in blocking:  # Prioritize moving blocking vehicles
                            priority -= 3
                        if vehicle == 'A' and direction > 0:  # Prioritize moving red car right
                            priority -= 2
                            
                        heapq.heappush(queue, (priority, new_board, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")