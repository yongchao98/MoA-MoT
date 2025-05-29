from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_target_state(board):
    # Check if white knights are at original black positions and vice versa
    return (board[0][1] == 'w' and board[0][3] == 'w' and 
            board[2][2] == 'B' and board[1][3] == 'B')

def find_solution():
    # Initial board state
    initial_board = [
        ['.', 'B', '.', 'B'],
        ['.', '.', '.', 'w'],
        ['.', '.', 'w', '.']
    ]
    
    # BFS queue: (board, moves_list, is_white_turn)
    queue = deque([(initial_board, [], True)])
    seen = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        piece = 'w' if is_white_turn else 'B'
        
        # Find all pieces of current color
        pieces = []
        for y in range(3):
            for x in range(4):
                if board[y][x] == piece:
                    pieces.append((x,y))
        
        # Try all possible moves for each piece
        for from_pos in pieces:
            fx, fy = from_pos
            for to_pos in get_knight_moves(from_pos):
                tx, ty = to_pos
                if board[ty][tx] == '.':  # If target square is empty
                    # Make move
                    new_board = copy.deepcopy(board)
                    new_board[fy][fx] = '.'
                    new_board[ty][tx] = piece
                    
                    board_str = board_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        move = f"{piece},{chr(65+fx)}{3-fy},{chr(65+tx)}{3-ty}"
                        new_moves = moves + [move]
                        
                        if is_target_state(new_board):
                            return new_moves
                        
                        queue.append((new_board, new_moves, not is_white_turn))
    
    return "No"

solution = find_solution()
print(solution)