from collections import deque

def parse_board(board_str):
    return [list(row) for row in board_str.strip().split('\n')]

def board_to_string(board):
    return '\n'.join(''.join(row) for row in board)

def get_vehicle_positions(board):
    vehicles = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] not in '.x':
                if board[i][j] not in vehicles:
                    vehicles[board[i][j]] = [(i, j)]
                else:
                    vehicles[board[i][j]].append((i, j))
    return vehicles

def is_horizontal(positions):
    return positions[0][0] == positions[1][0]

def make_move(board, vehicle, positions, direction):
    new_board = [row[:] for row in board]
    # Clear old positions
    for i, j in positions:
        new_board[i][j] = '.'
    # Set new positions
    for i, j in positions:
        if is_horizontal(positions):
            new_board[i][j + direction] = vehicle
        else:
            new_board[i + direction][j] = vehicle
    return new_board

def get_possible_moves(board, vehicles):
    moves = []
    for vehicle, positions in vehicles.items():
        horizontal = is_horizontal(positions)
        if horizontal:
            # Try moving left
            leftmost = min(p[1] for p in positions)
            if leftmost > 0 and board[positions[0][0]][leftmost-1] == '.':
                moves.append((vehicle, -1))
            # Try moving right
            rightmost = max(p[1] for p in positions)
            if rightmost < len(board[0])-1 and board[positions[0][0]][rightmost+1] == '.':
                moves.append((vehicle, 1))
        else:
            # Try moving up
            topmost = min(p[0] for p in positions)
            if topmost > 0 and board[topmost-1][positions[0][1]] == '.':
                moves.append((vehicle, -1))
            # Try moving down
            bottommost = max(p[0] for p in positions)
            if bottommost < len(board)-1 and board[bottommost+1][positions[0][1]] == '.':
                moves.append((vehicle, 1))
    return moves

def solve_puzzle():
    initial_board = parse_board("""..xBBx
...E..
AA.E..
.CCC..
..DD..
......""")
    
    queue = deque([(initial_board, [])])
    seen = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.popleft()
        
        # Check if solved (AA is at the exit)
        vehicles = get_vehicle_positions(current_board)
        if 'A' in vehicles and max(p[1] for p in vehicles['A']) == len(current_board[0])-1:
            return moves
        
        # Try all possible moves
        for vehicle, direction in get_possible_moves(current_board, vehicles):
            new_board = make_move(current_board, vehicle, vehicles[vehicle], direction)
            board_str = board_to_string(new_board)
            
            if board_str not in seen:
                seen.add(board_str)
                new_moves = moves + [f"{vehicle}{'+' if direction > 0 else '-'}{abs(direction)}"]
                queue.append((new_board, new_moves))
    
    return None

solution = solve_puzzle()
if solution:
    print("<<<" + " ".join(solution) + ">>>")
else:
    print("No solution found")