from collections import deque
import copy

def print_board(board):
    for i, row in enumerate(board):
        print(f"{3-i} |", end=" ")
        print(" | ".join(cell if cell != '.' else ' ' for cell in row), end=" ")
        print("|")
        print("  " + "-" * 19)
    print("    A   B   C   D")

def get_knight_moves(pos):
    row, col = pos
    moves = [
        (row+2, col+1), (row+2, col-1),
        (row-2, col+1), (row-2, col-1),
        (row+1, col+2), (row+1, col-2),
        (row-1, col+2), (row-1, col-2)
    ]
    return [(r,c) for r,c in moves if 0 <= r < 4 and 0 <= c < 4]

def algebraic_to_coords(pos):
    col = ord(pos[0]) - ord('A')
    row = 3 - (int(pos[1]) - 1)
    return (row, col)

def coords_to_algebraic(row, col):
    return f"{chr(col + ord('A'))}{3-row+1}"

def make_move(board, move):
    color, from_pos, to_pos = move.split(',')
    from_row, from_col = algebraic_to_coords(from_pos)
    to_row, to_col = algebraic_to_coords(to_pos)
    
    new_board = copy.deepcopy(board)
    new_board[from_row][from_col] = '.'
    new_board[to_row][to_col] = color
    return new_board

def is_valid_move(board, from_pos, to_pos, color):
    from_row, from_col = from_pos
    to_row, to_col = to_pos
    
    # Check if source has correct color piece
    if board[from_row][from_col] != color:
        return False
    
    # Check if destination is empty
    if board[to_row][to_col] != '.':
        return False
    
    # Check if move is L-shaped
    row_diff = abs(to_row - from_row)
    col_diff = abs(to_col - from_col)
    return (row_diff == 2 and col_diff == 1) or (row_diff == 1 and col_diff == 2)

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_puzzle():
    initial_board = [
        ['.', '.', 'B', '.'],
        ['.', 'B', '.', '.'],
        ['.', 'w', '.', 'w'],
        ['.', '.', '.', '.']
    ]
    
    print("Initial position:")
    print_board(initial_board)
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    seen = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        # Check if we've reached target positions
        if len(moves) > 20:  # Limit search depth
            continue
            
        piece = 'w' if is_white_turn else 'B'
        
        # Get all pieces of current color
        pieces = [(r,c) for r in range(4) for c in range(4) if board[r][c] == piece]
        
        for from_pos in pieces:
            valid_moves = get_knight_moves(from_pos)
            for to_pos in valid_moves:
                if board[to_pos[0]][to_pos[1]] == '.':
                    new_board = copy.deepcopy(board)
                    new_board[from_pos[0]][from_pos[1]] = '.'
                    new_board[to_pos[0]][to_pos[1]] = piece
                    
                    board_str = board_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        move = f"{piece},{coords_to_algebraic(*from_pos)},{coords_to_algebraic(*to_pos)}"
                        new_moves = moves + [move]
                        
                        # Check if we've reached the target state
                        if (piece == 'B' and 
                            set((r,c) for r in range(4) for c in range(4) if new_board[r][c] == 'w') == {(1,1), (0,2)} and
                            set((r,c) for r in range(4) for c in range(4) if new_board[r][c] == 'B') == {(2,1), (2,3)}):
                            return new_moves
                            
                        queue.append((new_board, new_moves, not is_white_turn))
    
    return "No"

solution = solve_puzzle()
print("\nSolution:", solution)