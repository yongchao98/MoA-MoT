from collections import deque

def get_knight_moves(pos):
    row, col = pos
    moves = [
        (row+2, col+1), (row+2, col-1),
        (row-2, col+1), (row-2, col-1),
        (row+1, col+2), (row+1, col-2),
        (row-1, col+2), (row-1, col-2)
    ]
    # Fix: Check bounds for 3x4 board
    return [(r, c) for r, c in moves if 0 <= r < 3 and 0 <= c < 4]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_target_reached(board):
    # White knights should be at A1 and D1, Black knights at A3 and B3
    return (board[0][0] == 'B' and board[0][1] == 'B' and 
            board[2][0] == 'w' and board[2][3] == 'w')

def find_solution():
    # Initial board: White knights at A3,B3; Black knights at A1,D1
    initial_board = [
        ['B', '.', '.', 'B'],  # Row 1
        ['.', '.', '.', '.'],  # Row 2
        ['w', 'w', '.', '.']   # Row 3
    ]
    
    queue = deque([(initial_board, [], True)])
    seen = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves, is_white_turn = queue.popleft()
        piece = 'w' if is_white_turn else 'B'
        
        # Find all pieces of current color
        pieces = []
        for r in range(3):
            for c in range(4):
                if current_board[r][c] == piece:
                    pieces.append((r, c))
        
        # Try all possible moves for each piece
        for from_pos in pieces:
            r, c = from_pos
            for to_r, to_c in get_knight_moves(from_pos):
                if current_board[to_r][to_c] == '.':
                    # Make move
                    new_board = [row[:] for row in current_board]
                    new_board[to_r][to_c] = piece
                    new_board[r][c] = '.'
                    
                    board_str = board_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        # Convert to chess notation (A1 = bottom-left)
                        from_notation = f"{chr(65+c)}{3-r}"
                        to_notation = f"{chr(65+to_c)}{3-to_r}"
                        move = f"{piece},{from_notation},{to_notation}"
                        new_moves = moves + [move]
                        
                        if is_target_reached(new_board):
                            return new_moves
                        
                        queue.append((new_board, new_moves, not is_white_turn))
    
    return None

solution = find_solution()
if solution:
    print(solution)
else:
    print("No")