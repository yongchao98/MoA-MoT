from collections import deque

def get_knight_moves(pos):
    row, col = pos
    moves = [
        (row+2, col+1), (row+2, col-1),
        (row-2, col+1), (row-2, col-1),
        (row+1, col+2), (row+1, col-2),
        (row-1, col+2), (row-1, col-2)
    ]
    return [(r, c) for r, c in moves if 0 <= r < 3 and 0 <= c < 4]

def pos_to_notation(pos):
    return f"{chr(65+pos[1])}{pos[0]+1}"

def get_pieces_positions(board):
    whites = []
    blacks = []
    empty = []
    for i in range(3):
        for j in range(4):
            if board[i][j] == 'w':
                whites.append((i, j))
            elif board[i][j] == 'B':
                blacks.append((i, j))
            elif board[i][j] == '.':
                empty.append((i, j))
    return whites, blacks, empty

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def find_solution():
    initial_board = [
        ['.', 'w', 'B', ' '],
        ['B', ' ', ' ', '.'],
        [' ', '.', 'w', ' ']
    ]
    
    target_positions = {
        'w': {(2, 1), (1, 0)},  # Original black positions
        'B': {(0, 1), (2, 2)}   # Original white positions
    }
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves, is_white_turn = queue.popleft()
        whites, blacks, empty = get_pieces_positions(current_board)
        
        # Check if we reached target state
        current_w_pos = set((r, c) for r, c in whites)
        current_b_pos = set((r, c) for r, c in blacks)
        if current_w_pos == target_positions['w'] and current_b_pos == target_positions['B']:
            return moves
        
        # Generate possible moves
        pieces = whites if is_white_turn else blacks
        piece_type = 'w' if is_white_turn else 'B'
        
        for piece in pieces:
            for new_pos in get_knight_moves(piece):
                if new_pos in [(r, c) for r, c in empty]:
                    # Make move
                    new_board = [row[:] for row in current_board]
                    new_board[piece[0]][piece[1]] = '.'
                    new_board[new_pos[0]][new_pos[1]] = piece_type
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        new_moves = moves + [f"{piece_type},{pos_to_notation(piece)},{pos_to_notation(new_pos)}"]
                        queue.append((new_board, new_moves, not is_white_turn))
    
    return None

solution = find_solution()
if solution:
    print(solution)
else:
    print("No")