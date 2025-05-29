from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 1 <= x <= 4 and 1 <= y <= 3]

def board_to_string(board):
    return ','.join(f"{pos}:{piece}" for pos, piece in sorted(board.items()))

def is_target_reached(board):
    # Check if white knights are at black's initial positions and vice versa
    white_positions = {pos for pos, piece in board.items() if piece == 'w'}
    black_positions = {pos for pos, piece in board.items() if piece == 'B'}
    return white_positions == {(3,3), (4,2)} and black_positions == {(1,1), (2,1)}

def solve_knight_swap():
    # Initial board state
    initial_board = {
        (1,1): 'w', (2,1): 'w',  # White knights
        (3,3): 'B', (4,2): 'B'   # Black knights
    }
    
    queue = deque([(initial_board, [], True)])  # (board, moves, is_white_turn)
    seen_states = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        if is_target_reached(board):
            return moves
            
        current_pieces = 'w' if is_white_turn else 'B'
        
        # Get all pieces of current color
        pieces = [(pos, piece) for pos, piece in board.items() if piece == current_pieces]
        
        # Try all possible moves for each piece
        for pos, piece in pieces:
            for new_pos in get_knight_moves(pos):
                if new_pos not in board:  # If target square is empty
                    # Make move
                    new_board = copy.deepcopy(board)
                    del new_board[pos]
                    new_board[new_pos] = piece
                    
                    board_state = board_to_string(new_board)
                    if board_state not in seen_states:
                        seen_states.add(board_state)
                        # Record move in chess notation (converting to A1 format)
                        from_square = f"{chr(pos[0]+64)}{pos[1]}"
                        to_square = f"{chr(new_pos[0]+64)}{new_pos[1]}"
                        new_moves = moves + [f"{piece},{from_square},{to_square}"]
                        queue.append((new_board, new_moves, not is_white_turn))
    
    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")