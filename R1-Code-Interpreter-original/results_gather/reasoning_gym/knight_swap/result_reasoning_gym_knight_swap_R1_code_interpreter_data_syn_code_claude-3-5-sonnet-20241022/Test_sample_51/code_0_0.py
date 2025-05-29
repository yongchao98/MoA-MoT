from collections import deque
import copy

def get_knight_moves(pos):
    x, y = ord(pos[0]) - ord('A'), int(pos[1]) - 1
    moves = [(2,1), (2,-1), (-2,1), (-2,-1), (1,2), (1,-2), (-1,2), (-1,-2)]
    valid = []
    for dx, dy in moves:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 4 and 0 <= new_y < 3:
            valid.append(chr(new_x + ord('A')) + str(new_y + 1))
    return valid

def board_to_string(board):
    return ','.join(f"{k}:{v}" for k, v in sorted(board.items()))

def is_goal_state(board):
    # Check if white knights are at original black positions and vice versa
    white_positions = {pos for pos, piece in board.items() if piece == 'w'}
    black_positions = {pos for pos, piece in board.items() if piece == 'B'}
    original_black = {'C1', 'D3'}
    original_white = {'A2', 'C2'}
    return white_positions == original_black and black_positions == original_white

def solve_knight_swap():
    # Initial board state
    initial_board = {
        'A2': 'w', 'C2': 'w',
        'C1': 'B', 'D3': 'B'
    }
    
    queue = deque([(initial_board, [], 'w')])  # board, moves, next_turn
    visited = {board_to_string(initial_board)}
    
    while queue:
        board, moves, turn = queue.popleft()
        
        if is_goal_state(board):
            return moves
        
        # Find pieces that can move
        pieces = [(pos, piece) for pos, piece in board.items() if piece == turn]
        
        for pos, piece in pieces:
            valid_moves = get_knight_moves(pos)
            for new_pos in valid_moves:
                if new_pos not in board:  # Empty square
                    new_board = copy.deepcopy(board)
                    del new_board[pos]
                    new_board[new_pos] = piece
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        new_moves = moves + [f"{piece},{pos},{new_pos}"]
                        queue.append((new_board, new_moves, 'B' if turn == 'w' else 'w'))
    
    return None

solution = solve_knight_swap()
if solution is None:
    print("No")
else:
    print(solution)