from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [
        (x+2, y+1), (x+2, y-1), (x-2, y+1), (x-2, y-1),
        (x+1, y+2), (x+1, y-2), (x-1, y+2), (x-1, y-2)
    ]
    return [(x, y) for x, y in moves if 1 <= x <= 4 and 1 <= y <= 3]

def board_to_string(board):
    return ','.join(f"{pos}:{piece}" for pos, piece in sorted(board.items()))

def is_target_state(board):
    target = {(3,3):'B', (3,1):'B', (3,2):'w', (3,4):'w'}
    return all(board.get((x,y)) == target.get((x,y)) for x,y in target)

def solve_knight_swap():
    # Initial board state
    initial_board = {
        (3,3): 'B', (3,4): 'w',
        (3,1): 'B', (3,2): '.',
        (2,2): '.', (2,1): '.'
    }
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves, is_white_turn = queue.popleft()
        
        if is_target_state(current_board):
            return moves
        
        # Find all pieces of current player
        current_pieces = [(pos, piece) for pos, piece in current_board.items() 
                         if piece == ('w' if is_white_turn else 'B')]
        
        # Try all possible moves for current pieces
        for pos, piece in current_pieces:
            for new_pos in get_knight_moves(pos):
                if current_board.get(new_pos) == '.':
                    # Make move
                    new_board = copy.deepcopy(current_board)
                    new_board[pos] = '.'
                    new_board[new_pos] = piece
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        move = f"{piece},{pos[1]}{pos[0]},{new_pos[1]}{new_pos[0]}"
                        queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

# Solve and print result
result = solve_knight_swap()
if result is None:
    print("<<<No>>>")
else:
    print(f"<<<{result}>>>")