from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_goal_state(board):
    # White knights should be at (A3, D3) - original black positions
    # Black knights should be at (A2, B3) - original white positions
    white_positions = {(0,0), (0,3)}  # A3, D3
    black_positions = {(1,0), (0,1)}  # A2, B3
    
    for y in range(3):
        for x in range(4):
            if board[y][x] == 'w' and (y,x) not in white_positions:
                return False
            if board[y][x] == 'B' and (y,x) not in black_positions:
                return False
    return True

def solve_knight_swap():
    initial_board = [
        ['B', 'w', '.', 'B'],
        ['w', '.', '.', '.'],
        ['.', '.', '.', '.']
    ]
    
    queue = deque([(initial_board, [], True)])
    visited = {board_to_string(initial_board)}
    max_moves = 20  # Prevent infinite loops
    
    while queue:
        board, moves, is_black_turn = queue.popleft()
        
        if len(moves) >= max_moves:
            continue
            
        pieces = []
        for i in range(3):
            for j in range(4):
                if board[i][j] == ('B' if is_black_turn else 'w'):
                    pieces.append((j, i))
        
        for piece in pieces:
            x, y = piece
            for new_x, new_y in get_knight_moves(piece):
                if board[new_y][new_x] == '.':
                    new_board = copy.deepcopy(board)
                    new_board[y][x] = '.'
                    new_board[new_y][new_x] = 'B' if is_black_turn else 'w'
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        move = f"{'B' if is_black_turn else 'w'},{chr(65+x)}{y+1},{chr(65+new_x)}{new_y+1}"
                        new_moves = moves + [move]
                        
                        if is_goal_state(new_board):
                            print(new_moves)
                            return
                        
                        queue.append((new_board, new_moves, not is_black_turn))
    
    print("No")

solve_knight_swap()