from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 4]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_goal_state(board):
    # Check if white knights are at B1,C1 and black knights at B2,C3
    return (board[0][2] == 'B' and board[1][1] == 'B' and 
            board[2][1] == 'w' and board[2][2] == 'w')

def get_piece_positions(board, piece):
    return [(i,j) for i in range(4) for j in range(4) if board[i][j] == piece]

def solve_knight_swap():
    # Initial board state (0,0 is top-left)
    initial_board = [
        ['.', '.', 'w', '.'],
        ['.', 'w', '.', '.'],
        ['.', 'B', 'B', '.']
    ]
    
    visited = set()
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        board_str = board_to_string(board)
        
        if board_str in visited:
            continue
        
        visited.add(board_str)
        
        if is_goal_state(board):
            return moves
        
        # Get current player's pieces
        piece = 'w' if is_white_turn else 'B'
        pieces = get_piece_positions(board, piece)
        
        # Try all possible moves for current player's pieces
        for from_pos in pieces:
            for to_pos in get_knight_moves(from_pos):
                i, j = to_pos
                if board[i][j] == '.':  # If target square is empty
                    # Make move
                    new_board = copy.deepcopy(board)
                    new_board[from_pos[0]][from_pos[1]] = '.'
                    new_board[i][j] = piece
                    
                    # Record move
                    from_square = f"{chr(from_pos[1]+65)}{3-from_pos[0]}"
                    to_square = f"{chr(j+65)}{3-i}"
                    move = f"{piece},{from_square},{to_square}"
                    
                    queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")