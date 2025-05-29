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
    # Check if white knights are at C3, C2 and black knights at B1, D2
    return (board[0][2] == 'w' and board[1][2] == 'w' and 
            board[2][1] == 'B' and board[1][3] == 'B')

def find_solution():
    # Initial board state
    initial_board = [
        ['.', '.', 'B', '.'],
        ['.', '.', 'B', 'w'],
        ['.', 'w', '.', '.']
    ]
    
    # Queue for BFS: (board, moves_list, is_black_turn)
    queue = deque([(initial_board, [], True)])
    visited = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_black_turn = queue.popleft()
        
        # Find current positions of pieces to move
        pieces = []
        for i in range(3):
            for j in range(4):
                if board[i][j] == ('B' if is_black_turn else 'w'):
                    pieces.append((j, i))  # x,y coordinates
        
        # Try all possible moves for current pieces
        for piece in pieces:
            for new_pos in get_knight_moves(piece):
                if board[new_pos[1]][new_pos[0]] == '.':
                    # Make move
                    new_board = copy.deepcopy(board)
                    new_board[piece[1]][piece[0]] = '.'
                    new_board[new_pos[1]][new_pos[0]] = 'B' if is_black_turn else 'w'
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        
                        # Record move
                        move = f"{'B' if is_black_turn else 'w'},{chr(65+piece[0])}{piece[1]+1},{chr(65+new_pos[0])}{new_pos[1]+1}"
                        new_moves = moves + [move]
                        
                        # Check if target reached
                        if is_target_state(new_board):
                            return new_moves
                        
                        queue.append((new_board, new_moves, not is_black_turn))
    
    return "No"

solution = find_solution()
print(solution)