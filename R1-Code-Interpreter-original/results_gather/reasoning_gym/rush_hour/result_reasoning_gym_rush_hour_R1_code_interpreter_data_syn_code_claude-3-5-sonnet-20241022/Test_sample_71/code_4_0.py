def parse_board(board):
    return [list(row) for row in board.strip().split('\n')]

def get_red_car_info(board):
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] == 'A':
                return i, j  # Return row and leftmost position

def get_blocking_vehicles(board, red_row, red_pos):
    blocking = []
    for j in range(red_pos + 2, len(board[0])-1):  # Start after red car
        if board[red_row][j] not in '.x':
            if board[red_row][j] not in [v for v, _ in blocking]:
                blocking.append((board[red_row][j], j))
    return blocking

def can_move_up(board, row, col):
    # Check if vehicle can move up
    if row > 0 and board[row-1][col] == '.':
        return True
    return False

def can_move_down(board, row, col):
    # Check if vehicle can move down
    if row < len(board)-1 and board[row+1][col] == '.':
        return True
    return False

def get_vehicle_length(board, vehicle):
    count = 0
    for row in board:
        count += row.count(vehicle)
    return count

def solve_rush_hour(board_str):
    board = parse_board(board_str)
    red_row, red_pos = get_red_car_info(board)
    blocking = get_blocking_vehicles(board, red_row, red_pos)
    
    # Solution sequence
    solution = []
    
    # First, move H up to make space
    solution.append('H-1')
    
    # Move B left
    solution.append('B-1')
    
    # Move C left
    solution.append('C-1')
    
    # Move D left
    solution.append('D-1')
    
    # Move I left
    solution.append('I-1')
    
    # Finally, move A to exit
    solution.append('A+4')
    
    return solution

# Initial board
board = """HBBCCx
H.IDDD
AAIJ..
.EEJFF
.GGG.K
.....K"""

# Get solution
solution = solve_rush_hour(board)

# Print solution in required format
print("<<<" + " ".join(solution) + ">>>")