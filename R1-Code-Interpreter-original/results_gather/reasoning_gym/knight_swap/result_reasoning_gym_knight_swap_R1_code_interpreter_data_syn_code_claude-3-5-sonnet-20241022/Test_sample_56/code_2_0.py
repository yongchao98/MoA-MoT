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

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_target_reached(board):
    # White knights should be at A1 and D1
    # Black knights should be at A3 and B3
    return (board[2][0] == 'B' and board[2][1] == 'B' and 
            board[0][0] == 'w' and board[0][3] == 'w')

def find_solution():
    initial_board = [
        ['B', '.', '.', 'B'],  # Row 1
        ['.', '.', '.', '.'],  # Row 2
        ['w', 'w', '.', '.']   # Row 3
    ]
    
    queue = deque([(initial_board, [], True)])
    seen = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves, is_white_turn = queue.popleft()
        
        if is_target_reached(current_board):
            return moves
            
        piece = 'w' if is_white_turn else 'B'
        pieces = [(r, c) for r in range(3) for c in range(4) 
                 if current_board[r][c] == piece]
        
        for r, c in pieces:
            for to_r, to_c in get_knight_moves((r, c)):
                if current_board[to_r][to_c] == '.':
                    new_board = [row[:] for row in current_board]
                    new_board[to_r][to_c] = piece
                    new_board[r][c] = '.'
                    
                    board_str = board_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        from_notation = f"{chr(65+c)}{3-r}"
                        to_notation = f"{chr(65+to_c)}{3-to_r}"
                        move = f"{piece},{from_notation},{to_notation}"
                        
                        queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

solution = find_solution()
print(solution if solution else "No")