from collections import deque

def parse_board(board_str):
    return [list(line) for line in board_str.splitlines()]

def get_vehicle_positions(board):
    vehicles = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] not in '.x':
                if board[i][j] not in vehicles:
                    vehicles[board[i][j]] = []
                vehicles[board[i][j]].append((i, j))
    
    # Determine orientation and sort positions
    for v in vehicles:
        pos = vehicles[v]
        vehicles[v] = {
            'positions': sorted(pos),
            'orientation': 'horizontal' if pos[0][0] == pos[1][0] else 'vertical'
        }
    return vehicles

def can_move(board, positions, di, dj):
    for i, j in positions:
        ni, nj = i + di, j + dj
        if not (0 <= ni < len(board) and 0 <= nj < len(board[0])):
            return False
        if board[ni][nj] not in ['.', board[i][j]]:
            return False
    return True

def move_vehicle(board, positions, di, dj):
    new_board = [row[:] for row in board]
    vehicle = board[positions[0][0]][positions[0][1]]
    
    for i, j in positions:
        new_board[i][j] = '.'
    
    for i, j in positions:
        new_board[i + di][j + dj] = vehicle
    
    return new_board

def board_to_string(board):
    return '\n'.join(''.join(row) for row in board)

def solve_rush_hour(board_str):
    board = parse_board(board_str)
    vehicles = get_vehicle_positions(board)
    
    queue = deque([(board, [])])
    seen = {board_to_string(board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if solved (red car at exit)
        red_car = [(i, j) for i in range(len(current_board)) 
                   for j in range(len(current_board[0])) if current_board[i][j] == 'A']
        if max(j for _, j in red_car) == len(current_board[0]) - 1:
            return moves
        
        # Try moves for each vehicle
        for vehicle, info in vehicles.items():
            current_pos = [(i, j) for i in range(len(current_board)) 
                          for j in range(len(current_board[0])) 
                          if current_board[i][j] == vehicle]
            
            if info['orientation'] == 'horizontal':
                # Try moving left/right
                for steps in [-1, 1]:
                    if can_move(current_board, current_pos, 0, steps):
                        new_board = move_vehicle(current_board, current_pos, 0, steps)
                        board_str = board_to_string(new_board)
                        if board_str not in seen:
                            seen.add(board_str)
                            new_moves = moves + [f"{vehicle}{'+' if steps > 0 else ''}{steps}"]
                            queue.append((new_board, new_moves))
            else:
                # Try moving up/down
                for steps in [-1, 1]:
                    if can_move(current_board, current_pos, steps, 0):
                        new_board = move_vehicle(current_board, current_pos, steps, 0)
                        board_str = board_to_string(new_board)
                        if board_str not in seen:
                            seen.add(board_str)
                            new_moves = moves + [f"{vehicle}{'+' if steps > 0 else ''}{steps}"]
                            queue.append((new_board, new_moves))
    
    return None

# Initial board
board = """BBBCCM
..JDDM
AAJK..
IEEKFF
IGGKLx
HH..L."""

solution = solve_rush_hour(board)
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")