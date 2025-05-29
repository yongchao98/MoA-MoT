from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 1 <= x <= 4 and 1 <= y <= 3]

def board_to_string(board):
    return ','.join(f"{pos}:{piece}" for pos,piece in sorted(board.items()))

def is_target_reached(board):
    # Check if white knights are at D2,D1 and black knights at A3,C1
    return (board.get((4,2)) == 'w' and board.get((4,1)) == 'w' and
            board.get((1,3)) == 'B' and board.get((3,1)) == 'B')

def find_solution():
    # Initial board state
    initial_board = {
        (1,3): 'w',  # A3
        (3,1): 'w',  # C1
        (4,2): 'B',  # D2
        (4,1): 'B',  # D1
    }
    
    # BFS queue: (board_state, move_history, is_white_turn)
    queue = deque([(initial_board, [], True)])
    seen = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        # Get all pieces of current player
        current_pieces = [(pos, piece) for pos, piece in board.items() 
                         if piece == ('w' if is_white_turn else 'B')]
        
        # Try all possible moves for current pieces
        for pos, piece in current_pieces:
            for new_pos in get_knight_moves(pos):
                if new_pos not in board:  # If target square is empty
                    # Make move
                    new_board = copy.deepcopy(board)
                    del new_board[pos]
                    new_board[new_pos] = piece
                    
                    board_str = board_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        
                        # Record move
                        move = f"{piece},{chr(pos[0]+64)}{pos[1]},{chr(new_pos[0]+64)}{new_pos[1]}"
                        new_moves = moves + [move]
                        
                        # Check if target reached
                        if is_target_reached(new_board):
                            return new_moves
                        
                        queue.append((new_board, new_moves, not is_white_turn))
    
    return None

# Find and print solution
solution = find_solution()
if solution:
    print(solution)
else:
    print("No")