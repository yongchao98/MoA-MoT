from heapq import heappush, heappop

def create_board():
    return [
        ['.', '.', 'x', 'B', 'B', 'B'],
        ['G', 'H', 'I', '.', 'C', 'C'],
        ['G', 'H', 'I', 'A', 'A', 'K'],
        ['D', 'D', 'D', 'J', '.', 'K'],
        ['.', '.', '.', 'J', 'E', 'E'],
        ['.', 'F', 'F', 'F', '.', 'x']
    ]

def get_blocking_cars(board):
    # Get cars blocking the red car's path
    blocking = set()
    red_row = 2
    red_col = None
    
    # Find red car position
    for j in range(6):
        if board[red_row][j] == 'A':
            red_col = j
            break
    
    # Find cars blocking the path
    for j in range(red_col + 2, 6):
        if board[red_row][j] not in '.x':
            blocking.add(board[red_row][j])
    
    return blocking

def can_move(board, row, col, car, is_horizontal, direction):
    if is_horizontal:
        new_col = col + direction
        return 0 <= new_col < 6 and board[row][new_col] == '.'
    else:
        new_row = row + direction
        return 0 <= new_row < 6 and board[new_row][col] == '.'

def make_move(board, move):
    new_board = [row[:] for row in board]
    car = move[0]
    direction = 1 if '+' in move else -1
    
    # Find car positions
    positions = []
    is_horizontal = None
    for i in range(6):
        for j in range(6):
            if new_board[i][j] == car:
                positions.append((i, j))
    
    positions.sort()
    is_horizontal = positions[0][0] == positions[1][0]
    
    # Move car
    if is_horizontal:
        row = positions[0][0]
        if direction > 0:  # Move right
            if can_move(board, row, positions[-1][1], car, True, 1):
                new_board[row][positions[-1][1] + 1] = car
                new_board[row][positions[0][1]] = '.'
        else:  # Move left
            if can_move(board, row, positions[0][1], car, True, -1):
                new_board[row][positions[0][1] - 1] = car
                new_board[row][positions[-1][1]] = '.'
    else:
        col = positions[0][1]
        if direction > 0:  # Move down
            if can_move(board, positions[-1][0], col, car, False, 1):
                new_board[positions[-1][0] + 1][col] = car
                new_board[positions[0][0]][col] = '.'
        else:  # Move up
            if can_move(board, positions[0][0], col, car, False, -1):
                new_board[positions[0][0] - 1][col] = car
                new_board[positions[-1][0]][col] = '.'
    
    return new_board

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_puzzle():
    initial_board = create_board()
    queue = [(initial_board, [])]
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves = queue.pop(0)
        blocking = get_blocking_cars(current_board)
        
        # If path is clear for red car
        if not blocking:
            # Find red car's rightmost position
            red_col = 0
            for j in range(6):
                if current_board[2][j] == 'A':
                    red_col = max(red_col, j)
            return moves + [f"A+{5 - red_col}"]
        
        # First try to move blocking cars
        for car in blocking:
            # Try moving up/down for vertical cars
            for direction in [-1, 1]:
                move = f"{car}{'+' if direction > 0 else '-'}1"
                new_board = make_move(current_board, move)
                board_str = board_to_string(new_board)
                if board_str != board_to_string(current_board) and board_str not in visited:
                    visited.add(board_str)
                    queue.append((new_board, moves + [move]))
        
        # Then try moving cars that might help clear the path
        for i in range(6):
            for j in range(6):
                if current_board[i][j] not in '.x':
                    car = current_board[i][j]
                    if car not in blocking and car != 'A':
                        for direction in [-1, 1]:
                            move = f"{car}{'+' if direction > 0 else '-'}1"
                            new_board = make_move(current_board, move)
                            board_str = board_to_string(new_board)
                            if board_str != board_to_string(current_board) and board_str not in visited:
                                visited.add(board_str)
                                queue.append((new_board, moves + [move]))

solution = solve_puzzle()
print('<<<' + ' '.join(solution) + '>>>')