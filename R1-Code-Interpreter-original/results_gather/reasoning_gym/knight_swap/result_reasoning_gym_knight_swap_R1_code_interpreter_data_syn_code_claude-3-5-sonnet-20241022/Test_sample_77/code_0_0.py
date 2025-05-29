from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def get_positions(board, piece):
    return [(i,j) for i in range(4) for j in range(3) 
            if board[j][i] == piece]

def make_move(board, start, end):
    new_board = [list(row) for row in board]
    piece = new_board[start[1]][start[0]]
    new_board[start[1]][start[0]] = '.'
    new_board[end[1]][end[0]] = piece
    return new_board

def solve_puzzle():
    initial_board = [
        ['B', '.', '.', ' '],
        [' ', ' ', 'B', 'w'],
        ['.', 'w', '.', '.']
    ]
    target_positions_white = get_positions([
        ['w', '.', '.', ' '],
        [' ', ' ', 'w', 'B'],
        ['.', 'B', '.', '.']
    ], 'w')
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    visited = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        current_pieces = get_positions(board, 'w' if is_white_turn else 'B')
        
        # Check if white pieces are in target positions when it's black's turn
        if not is_white_turn:
            white_positions = get_positions(board, 'w')
            if all(pos in target_positions_white for pos in white_positions):
                black_positions = get_positions(board, 'B')
                if all(pos not in target_positions_white for pos in black_positions):
                    return moves
        
        for start in current_pieces:
            for end in get_knight_moves(start):
                if board[end[1]][end[0]] == '.':
                    new_board = make_move(board, start, end)
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        move = f"{'w' if is_white_turn else 'B'},{chr(65+start[0])}{start[1]+1},{chr(65+end[0])}{end[1]+1}"
                        queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

solution = solve_puzzle()
if solution:
    print(solution)
else:
    print("No")