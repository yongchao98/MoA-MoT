def get_car_positions(board):
    cars = {}
    for i in range(len(board)):
        for j in range(len(board[i])):
            if board[i][j] != '.' and board[i][j] != 'x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = [(i, j)]
                else:
                    cars[board[i][j]].append((i, j))
    return cars

def try_move(board, car, direction):
    new_board = [list(row) for row in board]
    cars = get_car_positions(board)
    positions = cars[car]
    is_horizontal = positions[0][0] == positions[1][0]
    
    # Clear current positions
    for pos in positions:
        new_board[pos[0]][pos[1]] = '.'
    
    # Try new positions
    try:
        for pos in positions:
            new_row = pos[0] if is_horizontal else pos[0] + direction
            new_col = pos[1] + direction if is_horizontal else pos[1]
            if new_row < 0 or new_col < 0:
                return None
            if new_board[new_row][new_col] not in ['.']:
                return None
            new_board[new_row][new_col] = car
    except IndexError:
        return None
    
    return [''.join(row) for row in new_board]

def solve_direct(board):
    # Direct solution attempt based on the board layout
    moves = []
    current_board = board[:]
    
    # Step 1: Move K left
    new_board = try_move(current_board, 'K', -1)
    if new_board:
        moves.append('K-1')
        current_board = new_board
    else:
        return None
        
    # Step 2: Move D left
    new_board = try_move(current_board, 'D', -1)
    if new_board:
        moves.append('D-1')
        current_board = new_board
    else:
        return None
        
    # Step 3: Move J left
    new_board = try_move(current_board, 'J', -1)
    if new_board:
        moves.append('J-1')
        current_board = new_board
    else:
        return None
        
    # Step 4: Move A right
    new_board = try_move(current_board, 'A', 1)
    if new_board:
        moves.append('A+1')
        current_board = new_board
        return moves
    
    return None

initial_board = [
    "BBB.xK",
    "CCDDDK",
    "HAAJ..",
    "H..JEE",
    "..IFFx",
    "GGI..."
]

solution = solve_direct(initial_board)
if solution:
    print('<<<' + ' '.join(solution) + '>>>')
else:
    print("No solution found")