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
    # Check if white knights are at original black positions and vice versa
    return board[0][1] == '.' and board[1][3] == '.' and \
           board[2][2] == 'w' and board[2][3] == 'w' and \
           board[0][1] == '.' and board[1][3] == '.' and \
           'B' in [board[0][1], board[1][3]]

def find_solution():
    # Initial board state
    initial_board = [
        ['.', 'w', '.', '.'],
        ['.', '.', '.', 'w'],
        ['.', '.', 'B', 'B']
    ]
    
    # Queue for BFS: (board, moves_list, turn)
    queue = deque([(initial_board, [], True)])  # True for white's turn
    seen = {board_to_string(initial_board)}
    
    while queue:
        board, moves, white_turn = queue.popleft()
        
        # Find all pieces of current player
        pieces = []
        for i in range(3):
            for j in range(4):
                if board[i][j] == ('w' if white_turn else 'B'):
                    pieces.append((j, i))  # x,y coordinates
        
        # Try all possible moves for each piece
        for piece in pieces:
            for new_pos in get_knight_moves(piece):
                if board[new_pos[1]][new_pos[0]] == '.':  # If target square is empty
                    # Make move
                    new_board = copy.deepcopy(board)
                    new_board[piece[1]][piece[0]] = '.'
                    new_board[new_pos[1]][new_pos[0]] = 'w' if white_turn else 'B'
                    
                    board_str = board_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        
                        # Record move
                        move = f"{'w' if white_turn else 'B'},{chr(65+piece[0])}{piece[1]+1},{chr(65+new_pos[0])}{new_pos[1]+1}"
                        new_moves = moves + [move]
                        
                        # Check if target state reached
                        if new_board[2][2] == 'w' and new_board[2][3] == 'w' and \
                           new_board[0][1] == 'B' and new_board[1][3] == 'B':
                            return new_moves
                        
                        queue.append((new_board, new_moves, not white_turn))
    
    return "No"

# Solve and print result
result = find_solution()
print(result)