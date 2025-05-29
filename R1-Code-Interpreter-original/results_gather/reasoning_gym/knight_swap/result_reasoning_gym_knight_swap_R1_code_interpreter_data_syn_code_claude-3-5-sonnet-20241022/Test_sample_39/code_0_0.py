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
    # Check if black knights are at C2,D2 and white knights at B1,B2
    target_positions = {
        (2,1): 'w', (2,2): 'w',  # B1,B2
        (3,2): 'B', (4,2): 'B'   # C2,D2
    }
    return all(board.get((x,y)) == piece for (x,y),piece in target_positions.items())

def solve_knight_swap():
    # Initial board state
    initial_board = {
        (2,1): 'B', (2,2): 'B',  # B1,B2
        (3,2): 'w', (4,2): 'w'   # C2,D2
    }
    
    queue = deque([(initial_board, [], 'B')])  # (board, moves, next_player)
    seen = {board_to_string(initial_board)}
    
    while queue:
        board, moves, player = queue.popleft()
        
        if is_target_reached(board):
            return moves
            
        # Find all pieces of current player
        pieces = [(pos,piece) for pos,piece in board.items() if piece == player]
        
        # Try all possible moves for each piece
        for pos, piece in pieces:
            for new_pos in get_knight_moves(pos):
                if new_pos not in board:  # If target square is empty
                    # Make move
                    new_board = copy.deepcopy(board)
                    del new_board[pos]
                    new_board[new_pos] = piece
                    
                    board_str = board_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        new_moves = moves + [f"{piece},{pos[0]}{pos[1]},{new_pos[0]}{new_pos[1]}"]
                        queue.append((new_board, new_moves, 'w' if player == 'B' else 'B'))

    return None

# Solve and format the result
result = solve_knight_swap()
if result is None:
    print("<<<No>>>")
else:
    formatted_result = [move.replace(',', '') for move in result]
    print(f"<<<{formatted_result}>>>")