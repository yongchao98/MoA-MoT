from collections import deque
import copy

def get_knight_moves(pos):
    row, col = pos
    moves = [
        (row+2, col+1), (row+2, col-1),
        (row-2, col+1), (row-2, col-1),
        (row+1, col+2), (row+1, col-2),
        (row-1, col+2), (row-1, col-2)
    ]
    return [(r,c) for r,c in moves if 1 <= r <= 3 and 1 <= c <= 4]

def pos_to_notation(pos):
    return f"{chr(64+pos[1])}{pos[0]}"

def solve_knight_swap():
    # Initial board state (1-based indexing for rows and columns)
    initial_board = {
        (3,1): 'B', (3,2): '.', (3,3): 'w', (3,4): '.',
        (2,1): '.', (2,2): '.', (2,3): '.', (2,4): '.',
        (1,1): '.', (1,2): 'w', (1,3): '.', (1,4): 'B'
    }
    
    # Target positions
    white_targets = {(3,1), (1,4)}  # Where black knights start
    black_targets = {(3,3), (1,2)}  # Where white knights start
    
    # BFS queue: (board_state, moves_list, is_white_turn)
    queue = deque([(initial_board, [], True)])
    seen = set()
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        # Check if we've reached target state
        whites = {pos for pos, piece in board.items() if piece == 'w'}
        blacks = {pos for pos, piece in board.items() if piece == 'B'}
        if whites == white_targets and blacks == black_targets:
            return moves
        
        # Get current positions of pieces to move
        current_pieces = whites if is_white_turn else blacks
        piece_char = 'w' if is_white_turn else 'B'
        
        # Try all possible moves for current player
        for pos in current_pieces:
            for new_pos in get_knight_moves(pos):
                if new_pos in board and board[new_pos] == '.':  # Check if move is valid
                    # Make move
                    new_board = copy.deepcopy(board)
                    new_board[pos] = '.'
                    new_board[new_pos] = piece_char
                    
                    # Create board hash for seen positions
                    board_hash = str(new_board) + str(not is_white_turn)
                    if board_hash in seen:
                        continue
                    seen.add(board_hash)
                    
                    # Add move
                    move = f"{piece_char},{pos_to_notation(pos)},{pos_to_notation(new_pos)}"
                    new_moves = moves + [move]
                    
                    queue.append((new_board, new_moves, not is_white_turn))
    
    return None

solution = solve_knight_swap()
if solution is None:
    print("<<<No>>>")
else:
    print(f"<<<{solution}>>>")