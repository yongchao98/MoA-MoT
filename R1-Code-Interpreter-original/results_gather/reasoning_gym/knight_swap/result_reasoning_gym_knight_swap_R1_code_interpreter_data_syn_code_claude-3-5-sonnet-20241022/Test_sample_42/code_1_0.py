from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def get_piece_positions(board, piece):
    return [(x,y) for y in range(3) for x in range(4) if board[y][x] == piece]

def is_target_reached(board):
    white_pos = set(get_piece_positions(board, 'w'))
    black_pos = set(get_piece_positions(board, 'B'))
    target_white = {(0,0), (3,1)}  # A3, D2 positions for white
    target_black = {(1,0), (2,0)}  # B3, C3 positions for black
    return white_pos == target_white and black_pos == target_black

def make_move(board, start, end):
    new_board = [list(row) for row in board]
    piece = new_board[start[1]][start[0]]
    new_board[start[1]][start[0]] = '.'
    new_board[end[1]][end[0]] = piece
    return new_board

def solve_puzzle():
    # Corrected initial board (0,0 is top-left)
    initial_board = [
        ['B','w','w','.'],
        ['.','.','.','B'],
        ['.','.','.','.']
    ]
    
    visited = set()
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        board_str = board_to_string(board)
        
        if board_str in visited:
            continue
        visited.add(board_str)
        
        if is_target_reached(board):
            return moves
            
        current_pieces = get_piece_positions(board, 'w' if is_white_turn else 'B')
        for start in current_pieces:
            for end in get_knight_moves(start):
                if board[end[1]][end[0]] == '.':
                    new_board = make_move(board, start, end)
                    # Convert coordinates to chess notation (A3 = (0,0))
                    start_pos = f"{chr(65+start[0])}{3-start[1]}"
                    end_pos = f"{chr(65+end[0])}{3-end[1]}"
                    move = f"{'w' if is_white_turn else 'B'},{start_pos},{end_pos}"
                    queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

solution = solve_puzzle()
print("<<<" + str(solution) + ">>>" if solution else "<<<No>>>")