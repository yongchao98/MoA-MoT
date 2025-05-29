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
    return [(r,c) for r,c in moves if 0 <= r < 3 and 0 <= c < 3]

def pos_to_notation(pos):
    return chr(pos[1] + ord('A')) + str(3 - pos[0])

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_knight_swap():
    # Initial board
    initial_board = [
        ['.',  'w', '.'],
        ['.', '.', 'B'],
        ['B', '.', 'w']
    ]
    
    # Target board
    target_board = [
        ['.', 'B', '.'],
        ['.', '.', 'w'],
        ['w', '.', 'B']
    ]
    
    queue = deque([(initial_board, [], True)])  # board, moves, white_turn
    seen = {board_to_string(initial_board)}
    
    while queue:
        board, moves, white_turn = queue.popleft()
        
        # Check if we reached target
        if board == target_board:
            return moves
            
        # Find all pieces that can move
        piece = 'w' if white_turn else 'B'
        for r in range(3):
            for c in range(3):
                if board[r][c] == piece:
                    # Try all possible moves for this piece
                    for new_r, new_c in get_knight_moves((r,c)):
                        if board[new_r][new_c] == '.':
                            # Make move
                            new_board = copy.deepcopy(board)
                            new_board[r][c] = '.'
                            new_board[new_r][new_c] = piece
                            
                            board_str = board_to_string(new_board)
                            if board_str not in seen:
                                seen.add(board_str)
                                move = f"{piece},{pos_to_notation((r,c))},{pos_to_notation((new_r,new_c))}"
                                queue.append((new_board, moves + [move], not white_turn))
    
    return None

solution = solve_knight_swap()
if solution is None:
    print("<<<No>>>")
else:
    print(f"<<<{solution}>>>")