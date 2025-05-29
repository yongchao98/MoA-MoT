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
    return f"{chr(65+pos[1])}{3-pos[0]}"

def print_board(board):
    for row in board:
        print(row)
    print()

def find_solution():
    # Initial board state (0,0 is top-left)
    initial_board = [
        ['.', 'w', 'B', ' '],  # row 3
        ['B', ' ', ' ', '.'],  # row 2
        [' ', '.', 'w', ' ']   # row 1
    ]
    
    # Define target positions for each color
    target_w_positions = {(0, 2), (1, 0)}  # Where black knights started
    target_b_positions = {(0, 1), (2, 2)}  # Where white knights started
    
    queue = deque([(initial_board, [], True)])  # (board, moves_list, is_white_turn)
    visited = {str(initial_board)}
    
    while queue:
        current_board, moves, is_white_turn = queue.popleft()
        
        # Get current positions
        whites = []
        blacks = []
        empty = []
        for i in range(3):
            for j in range(4):
                if current_board[i][j] == 'w':
                    whites.append((i, j))
                elif current_board[i][j] == 'B':
                    blacks.append((i, j))
                elif current_board[i][j] == '.':
                    empty.append((i, j))
        
        # Check if we reached target state
        current_w_positions = set(whites)
        current_b_positions = set(blacks)
        
        if current_w_positions == target_w_positions and current_b_positions == target_b_positions:
            return moves
        
        # Generate moves for current player
        current_pieces = whites if is_white_turn else blacks
        piece_symbol = 'w' if is_white_turn else 'B'
        
        for piece_pos in current_pieces:
            valid_moves = get_knight_moves(piece_pos)
            for new_pos in valid_moves:
                if new_pos in empty:
                    # Make the move
                    new_board = [row[:] for row in current_board]
                    new_board[piece_pos[0]][piece_pos[1]] = '.'
                    new_board[new_pos[0]][new_pos[1]] = piece_symbol
                    
                    board_str = str(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        new_moves = moves + [f"{piece_symbol},{pos_to_notation(piece_pos)},{pos_to_notation(new_pos)}"]
                        queue.append((new_board, new_moves, not is_white_turn))
    
    return None

# Run the solution
solution = find_solution()
if solution:
    print(solution)
else:
    print("No")