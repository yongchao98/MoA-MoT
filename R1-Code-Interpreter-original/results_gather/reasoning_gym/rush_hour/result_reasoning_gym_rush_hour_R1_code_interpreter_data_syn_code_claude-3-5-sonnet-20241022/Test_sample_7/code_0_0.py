from collections import deque
import copy

def print_board(board):
    for row in board:
        print(''.join(row))

def is_solved(board):
    # Check if AA car can reach exit
    row = board[2]  # AA is in row 3 (index 2)
    pos = row.index('A')
    return all(c in '.A' for c in row[pos:])

def get_moves(board):
    moves = []
    cars = {}
    # Find all cars and their orientations
    for i in range(6):
        for j in range(6):
            if board[i][j].isalpha() and board[i][j] != 'x':
                if board[i][j] not in cars:
                    cars[board[i][j]] = [(i,j)]
                else:
                    cars[board[i][j]].append((i,j))
    
    # Try moves for each car
    for car in cars:
        coords = cars[car]
        if len(coords) >= 2:
            # Check if vertical
            if coords[0][1] == coords[1][1]:
                # Try moving up
                top = min(c[0] for c in coords)
                if top > 0 and board[top-1][coords[0][1]] == '.':
                    moves.append((car, -1))
                # Try moving down
                bottom = max(c[0] for c in coords)
                if bottom < 5 and board[bottom+1][coords[0][1]] == '.':
                    moves.append((car, 1))
            # Check if horizontal
            else:
                # Try moving left
                left = min(c[1] for c in coords)
                if left > 0 and board[coords[0][0]][left-1] == '.':
                    moves.append((car, -1))
                # Try moving right
                right = max(c[1] for c in coords)
                if right < 5 and board[coords[0][0]][right+1] == '.':
                    moves.append((car, 1))
    return moves

def apply_move(board, move):
    car, direction = move
    new_board = [list(row) for row in board]
    # Find car coordinates
    coords = []
    for i in range(6):
        for j in range(6):
            if new_board[i][j] == car:
                coords.append((i,j))
    
    # Check if vertical or horizontal
    if coords[0][1] == coords[1][1]:  # vertical
        coords.sort(key=lambda x: x[0])
        if direction < 0:  # move up
            new_board[coords[0][0]-1][coords[0][1]] = car
            new_board[coords[-1][0]][coords[-1][1]] = '.'
        else:  # move down
            new_board[coords[-1][0]+1][coords[-1][1]] = car
            new_board[coords[0][0]][coords[0][1]] = '.'
    else:  # horizontal
        coords.sort(key=lambda x: x[1])
        if direction < 0:  # move left
            new_board[coords[0][0]][coords[0][1]-1] = car
            new_board[coords[-1][0]][coords[-1][1]] = '.'
        else:  # move right
            new_board[coords[-1][0]][coords[-1][1]+1] = car
            new_board[coords[0][0]][coords[0][1]] = '.'
    
    return [''.join(row) for row in new_board]

def solve_puzzle():
    initial_board = [
        'E..x..',
        'EBBG..',
        'AAFG..',
        '..FCCC',
        '....H.',
        'xDDDH.'
    ]
    
    visited = set()
    queue = deque([(initial_board, [])])
    visited.add(tuple(initial_board))
    
    while queue:
        current_board, moves = queue.popleft()
        
        if is_solved(current_board):
            return ' '.join(f"{car}{'+' if dir > 0 else ''}{dir}" for car, dir in moves)
        
        for move in get_moves(current_board):
            new_board = apply_move(current_board, move)
            board_tuple = tuple(new_board)
            
            if board_tuple not in visited:
                visited.add(board_tuple)
                queue.append((new_board, moves + [move]))
    
    return "No solution found"

print(solve_puzzle())