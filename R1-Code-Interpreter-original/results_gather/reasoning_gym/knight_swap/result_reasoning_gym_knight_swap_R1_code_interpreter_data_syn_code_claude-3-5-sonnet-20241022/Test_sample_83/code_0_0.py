from collections import deque
import copy

def get_knight_moves(pos):
    x, y = ord(pos[0]) - ord('A'), int(pos[1]) - 1
    moves = []
    for dx, dy in [(2,1), (2,-1), (-2,1), (-2,-1), (1,2), (1,-2), (-1,2), (-1,-2)]:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 4 and 0 <= new_y < 3:
            moves.append(chr(new_x + ord('A')) + str(new_y + 1))
    return moves

def board_to_string(board):
    return ','.join(f"{pos}:{piece}" for pos, piece in sorted(board.items()))

def is_target_reached(board):
    # Check if white knights are at black's initial positions and vice versa
    white_positions = {pos for pos, piece in board.items() if piece == 'w'}
    black_positions = {pos for pos, piece in board.items() if piece == 'B'}
    return white_positions == {'C1', 'D2'} and black_positions == {'B3', 'C3'}

def find_solution():
    # Initial board setup
    initial_board = {
        'B3': 'w', 'C3': 'w',
        'C1': 'B', 'D2': 'B'
    }
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    visited = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        if is_target_reached(board):
            return moves
            
        current_pieces = 'w' if is_white_turn else 'B'
        current_positions = [pos for pos, piece in board.items() if piece == current_pieces]
        
        for pos in current_positions:
            for new_pos in get_knight_moves(pos):
                if new_pos not in board:  # Empty square
                    new_board = copy.deepcopy(board)
                    new_board[new_pos] = new_board.pop(pos)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in visited:
                        visited.add(board_str)
                        new_moves = moves + [f"{current_pieces},{pos},{new_pos}"]
                        queue.append((new_board, new_moves, not is_white_turn))
    
    return None

solution = find_solution()
if solution:
    print(solution)
else:
    print("No")