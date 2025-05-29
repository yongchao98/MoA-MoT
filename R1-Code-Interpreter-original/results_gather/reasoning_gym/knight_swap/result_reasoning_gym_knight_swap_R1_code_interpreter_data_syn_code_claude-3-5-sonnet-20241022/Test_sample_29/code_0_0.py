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

def is_goal_state(board):
    goal_whites = {(2,2), (4,2)}  # B2, D2
    goal_blacks = {(2,3), (2,2)}  # B3, B2
    current_whites = {pos for pos, piece in board.items() if piece == 'w'}
    current_blacks = {pos for pos, piece in board.items() if piece == 'B'}
    return current_whites == goal_whites and current_blacks == goal_blacks

def solve_knight_swap():
    # Initial board state
    initial_board = {
        (2,3): 'w',  # B3
        (2,2): 'w',  # B2
        (2,2): 'B',  # B2
        (4,2): 'B',  # D2
    }
    
    visited = set()
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        board_state = board_to_string(board)
        
        if board_state in visited:
            continue
        
        visited.add(board_state)
        
        if is_goal_state(board):
            return moves
        
        current_pieces = 'w' if is_white_turn else 'B'
        current_positions = [(pos, piece) for pos, piece in board.items() if piece == current_pieces]
        
        for pos, piece in current_positions:
            valid_moves = get_knight_moves(pos)
            for new_pos in valid_moves:
                if new_pos not in board:  # Empty square
                    new_board = copy.deepcopy(board)
                    del new_board[pos]
                    new_board[new_pos] = piece
                    move = f"{piece},{chr(pos[0]+64)}{pos[1]},{chr(new_pos[0]+64)}{new_pos[1]}"
                    queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

solution = solve_knight_swap()
if solution is None:
    print("<<<No>>>")
else:
    print(f"<<<{solution}>>>")