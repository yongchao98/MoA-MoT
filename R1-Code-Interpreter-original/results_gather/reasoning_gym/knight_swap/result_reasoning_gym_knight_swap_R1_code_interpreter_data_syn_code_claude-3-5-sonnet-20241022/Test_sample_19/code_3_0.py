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

def print_board(board):
    for row in board:
        print(row)
    print()

def is_target_state(board):
    # Check if white knights are at C1, D1 and black knights are at B3, D2
    return (board[2][2] == 'w' and board[2][3] == 'w' and  # C1, D1
            board[0][1] == 'B' and board[1][3] == 'B')     # B3, D2

def find_solution():
    # Initial board state
    initial_board = [
        ['.', 'w', '.', '.'],  # row 3 (index 0)
        ['.', '.', '.', 'w'],  # row 2 (index 1)
        ['.', '.', 'B', 'B']   # row 1 (index 2)
    ]
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    seen = {board_to_string(initial_board)}
    
    while queue:
        board, moves, white_turn = queue.popleft()
        
        if len(moves) > 10:  # Limit search depth
            continue
            
        # Find current player's pieces
        pieces = []
        for i in range(3):
            for j in range(4):
                if board[i][j] == ('w' if white_turn else 'B'):
                    pieces.append((j, i))  # (x, y)
        
        # Try all possible moves for each piece
        for piece in pieces:
            for new_pos in get_knight_moves(piece):
                if board[new_pos[1]][new_pos[0]] == '.':  # If target square is empty
                    # Make move
                    new_board = copy.deepcopy(board)
                    new_board[piece[1]][piece[0]] = '.'
                    new_board[new_pos[1]][new_pos[0]] = 'w' if white_turn else 'B'
                    
                    # Convert coordinates to chess notation
                    from_sq = f"{chr(65+piece[0])}{3-piece[1]}"
                    to_sq = f"{chr(65+new_pos[0])}{3-new_pos[1]}"
                    move = f"{'w' if white_turn else 'B'},{from_sq},{to_sq}"
                    
                    board_str = board_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        new_moves = moves + [move]
                        
                        # Check if we reached target state
                        if (white_turn and 
                            new_board[2][2] == 'w' and new_board[2][3] == 'w' and
                            new_board[0][1] == 'B' and new_board[1][3] == 'B'):
                            return new_moves
                        
                        queue.append((new_board, new_moves, not white_turn))
    
    return "No"

result = find_solution()
print(result)