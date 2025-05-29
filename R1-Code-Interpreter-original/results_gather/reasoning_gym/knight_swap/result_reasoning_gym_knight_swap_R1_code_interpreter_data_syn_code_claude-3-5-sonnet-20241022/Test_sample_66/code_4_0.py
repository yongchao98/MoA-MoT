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
    return f"{chr(65+pos[1])}{3-pos[0]}"

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def find_solution():
    initial_board = [
        ['.', '.', 'w'],  # row 3
        ['w', '.', '.'],  # row 2
        ['B', '.', 'B']   # row 1
    ]
    
    queue = deque([(initial_board, [], True)])
    seen = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        # Check if we've reached a winning position
        if (board[2][0] == 'w' and board[2][2] == 'w' and 
            board[0][2] == 'B' and board[1][0] == 'B'):
            return moves
        
        # Find current player's pieces
        piece = 'w' if is_white_turn else 'B'
        pieces = []
        for i in range(3):
            for j in range(3):
                if board[i][j] == piece:
                    pieces.append((i,j))
        
        # Try all possible moves
        for start_pos in pieces:
            for end_pos in get_knight_moves(start_pos):
                if board[end_pos[0]][end_pos[1]] == '.':
                    # Create new board state
                    new_board = [row[:] for row in board]
                    new_board[end_pos[0]][end_pos[1]] = piece
                    new_board[start_pos[0]][start_pos[1]] = '.'
                    
                    board_str = board_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        move = f"{piece},{pos_to_notation(start_pos)},{pos_to_notation(end_pos)}"
                        
                        # Prioritize moves that advance pieces towards their goals
                        if (piece == 'w' and end_pos[0] > start_pos[0]) or \
                           (piece == 'B' and end_pos[0] < start_pos[0]):
                            queue.appendleft((new_board, moves + [move], not is_white_turn))
                        else:
                            queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

solution = find_solution()
if solution is None:
    print("<<<No>>>")
else:
    print(f"<<<{solution}>>>")