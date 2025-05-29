from collections import deque
import copy

def get_knight_moves(pos):
    moves = [(2,1), (2,-1), (-2,1), (-2,-1), (1,2), (1,-2), (-1,2), (-1,-2)]
    x, y = pos
    valid = []
    for dx, dy in moves:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 4 and 0 <= new_y < 4:  # 4x4 board
            valid.append((new_x, new_y))
    return valid

def pos_to_notation(pos):
    return chr(pos[1] + ord('A')) + str(pos[0] + 1)

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def find_solution():
    # Initial board
    initial = [
        ['B', '.', '.', '.'],
        ['B', '.', '.', '.'],
        ['.', '.', 'w', 'w']
    ]
    
    # Target board
    target = [
        ['.', '.', '.', '.'],
        ['.', '.', '.', '.'],
        ['B', 'B', 'w', 'w']
    ]
    
    # Queue for BFS: (board, moves_list, is_white_turn)
    queue = deque([(initial, [], True)])
    visited = {board_to_string(initial)}
    
    while queue:
        board, moves, white_turn = queue.popleft()
        
        # Find current positions of all knights
        knights = []
        for i in range(3):
            for j in range(4):
                if board[i][j] in ['w', 'B']:
                    knights.append((i, j, board[i][j]))
        
        # Try all possible moves for current player
        for i, j, piece in knights:
            if (piece == 'w') == white_turn:  # Current player's pieces
                for new_i, new_j in get_knight_moves((i, j)):
                    if board[new_i][new_j] == '.':  # Empty square
                        # Make move
                        new_board = copy.deepcopy(board)
                        new_board[i][j] = '.'
                        new_board[new_i][new_j] = piece
                        
                        board_str = board_to_string(new_board)
                        if board_str not in visited:
                            visited.add(board_str)
                            
                            # Record move
                            move = f"{piece},{pos_to_notation((i,j))},{pos_to_notation((new_i,new_j))}"
                            new_moves = moves + [move]
                            
                            # Check if target reached
                            if new_board[2][0] == 'B' and new_board[2][1] == 'B' and \
                               new_board[2][2] == 'w' and new_board[2][3] == 'w':
                                return new_moves
                            
                            queue.append((new_board, new_moves, not white_turn))
    
    return "No"

solution = find_solution()
print(solution)