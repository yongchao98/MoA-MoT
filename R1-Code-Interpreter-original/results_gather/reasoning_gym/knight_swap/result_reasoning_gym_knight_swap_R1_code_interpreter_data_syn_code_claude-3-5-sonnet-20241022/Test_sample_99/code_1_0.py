from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [
        (x+2, y+1), (x+2, y-1), (x-2, y+1), (x-2, y-1),
        (x+1, y+2), (x+1, y-2), (x-1, y+2), (x-1, y-2)
    ]
    return [(x,y) for x,y in moves if 0 <= x < 3 and 0 <= y < 3]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def get_positions(board, piece):
    return [(i,j) for i in range(3) for j in range(3) if board[i][j] == piece]

def is_target_state(board):
    # White knights should be at B2, B1 (1,1 and 2,1)
    # Black knights should be at A3, A1 (0,0 and 2,0)
    white_pos = sorted(get_positions(board, 'w'))
    black_pos = sorted(get_positions(board, 'B'))
    return (white_pos == [(1,1), (2,1)] and 
            black_pos == [(0,0), (2,0)])

def solve_knight_swap():
    # Initialize board (0,0 is top-left)
    initial_board = [
        ['w', '.', '.'],  # A3,B3,C3
        ['B', '.', '.'],  # A2,B2,C2
        ['w', 'B', '.']   # A1,B1,C1
    ]
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    visited = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        if is_target_state(board):
            return moves
            
        piece = 'w' if is_white_turn else 'B'
        positions = get_positions(board, piece)
        
        for old_i, old_j in positions:
            for new_i, new_j in get_knight_moves((old_i, old_j)):
                if board[new_i][new_j] == '.':
                    new_board = copy.deepcopy(board)
                    new_board[new_i][new_j] = piece
                    new_board[old_i][old_j] = '.'
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        # Convert to chess notation (A1 is bottom-left)
                        from_sq = f"{chr(old_j+65)}{3-old_i}"
                        to_sq = f"{chr(new_j+65)}{3-new_i}"
                        move = f"{piece},{from_sq},{to_sq}"
                        queue.append((new_board, moves + [move], not is_white_turn))
    
    return "No"

print(solve_knight_swap())