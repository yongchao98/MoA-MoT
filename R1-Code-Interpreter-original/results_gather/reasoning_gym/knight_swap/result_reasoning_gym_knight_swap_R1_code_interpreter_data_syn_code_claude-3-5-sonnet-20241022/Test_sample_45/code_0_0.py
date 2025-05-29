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
    return [(r,c) for r,c in moves if 0 <= r < 4 and 0 <= c < 4]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_target_reached(board):
    # Check if white knights are at original black positions and vice versa
    white_positions = {(r,c) for r in range(4) for c in range(4) if board[r][c] == 'w'}
    target_white = {(2,1), (2,2)}  # B2, C3
    black_positions = {(r,c) for r in range(4) for c in range(4) if board[r][c] == 'B'}
    target_black = {(0,1), (0,3)}  # B1, D1
    return white_positions == target_white and black_positions == target_black

def get_positions(board, piece):
    return [(r,c) for r in range(4) for c in range(4) if board[r][c] == piece]

def solve_knight_swap():
    # Initial board
    initial_board = [
        ['.', '.', 'B', '.'],
        ['.', 'B', '.', '.'],
        ['.', 'w', '.', 'w'],
        ['.', '.', '.', '.']
    ]
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    seen = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        if is_target_reached(board):
            return moves
        
        piece = 'w' if is_white_turn else 'B'
        positions = get_positions(board, piece)
        
        for pos in positions:
            for new_pos in get_knight_moves(pos):
                r, c = new_pos
                if board[r][c] == '.':  # Can move to empty square
                    new_board = copy.deepcopy(board)
                    # Make move
                    new_board[pos[0]][pos[1]] = '.'
                    new_board[r][c] = piece
                    
                    board_str = board_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        # Convert position to algebraic notation
                        from_sq = f"{chr(pos[1]+65)}{3-pos[0]}"
                        to_sq = f"{chr(c+65)}{3-r}"
                        move = f"{piece},{from_sq},{to_sq}"
                        queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

# Solve and print result
solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")