from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [
        (x+2, y+1), (x+2, y-1), (x-2, y+1), (x-2, y-1),
        (x+1, y+2), (x+1, y-2), (x-1, y+2), (x-1, y-2)
    ]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_target_state(board):
    # Check if white knights are at target positions (where black knights started)
    white_positions = [(i,j) for i in range(4) for j in range(3) 
                      if board[j][i] == 'w']
    target_positions = [(2,2), (3,1)]  # C3, D2 in 0-based indexing
    return sorted(white_positions) == sorted(target_positions)

def solve_knight_swap():
    # Correct initial board state (0-based indexing)
    initial_board = [
        ['.', '.', 'w', 'B'],  # row 3
        ['.', '.', '.', 'B'],  # row 2
        ['w', '.', '.', '.']   # row 1
    ]
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_black_turn
    seen = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_black_turn = queue.popleft()
        
        # Find all pieces that can move
        pieces = []
        for i in range(4):
            for j in range(3):
                if board[j][i] == ('B' if is_black_turn else 'w'):
                    pieces.append((i,j))
        
        # Try all possible moves for each piece
        for piece in pieces:
            for new_pos in get_knight_moves(piece):
                if board[new_pos[1]][new_pos[0]] == '.':
                    # Make move
                    new_board = copy.deepcopy(board)
                    piece_type = new_board[piece[1]][piece[0]]
                    new_board[piece[1]][piece[0]] = '.'
                    new_board[new_pos[1]][new_pos[0]] = piece_type
                    
                    board_str = board_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        from_sq = f"{chr(piece[0]+65)}{3-piece[1]}"
                        to_sq = f"{chr(new_pos[0]+65)}{3-new_pos[1]}"
                        new_moves = moves + [f"{'B' if is_black_turn else 'w'},{from_sq},{to_sq}"]
                        
                        if is_target_state(new_board):
                            return new_moves
                        
                        queue.append((new_board, new_moves, not is_black_turn))
    
    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")