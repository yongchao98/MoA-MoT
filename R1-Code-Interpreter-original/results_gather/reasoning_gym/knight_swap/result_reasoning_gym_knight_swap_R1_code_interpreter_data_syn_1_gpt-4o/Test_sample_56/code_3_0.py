from collections import deque

def is_valid_move(x, y, board):
    return 0 <= x < 4 and 0 <= y < 3 and board[y][x] == '.'

def knight_moves(x, y):
    return [(x+2, y+1), (x+2, y-1), (x-2, y+1), (x-2, y-1),
            (x+1, y+2), (x+1, y-2), (x-1, y+2), (x-1, y-2)]

def board_to_tuple(board):
    return tuple(tuple(row) for row in board)

def swap_knights_bfs(initial_board):
    target_board = [['w', '.', '.', 'w'], ['.', '.', '.', '.'], ['B', 'B', '.', '.']]
    target_tuple = board_to_tuple(target_board)
    
    queue = deque([(initial_board, [], 'w')])
    visited = set()
    visited.add(board_to_tuple(initial_board))
    
    while queue:
        board, moves, turn = queue.popleft()
        
        if board_to_tuple(board) == target_tuple:
            return moves
        
        # Find current knight positions
        knight_positions = {'w': [], 'B': []}
        for y in range(3):
            for x in range(4):
                if board[y][x] in knight_positions:
                    knight_positions[board[y][x]].append((x, y))
        
        current_knights = knight_positions[turn]
        
        for i, (x, y) in enumerate(current_knights):
            for nx, ny in knight_moves(x, y):
                if is_valid_move(nx, ny, board):
                    # Make the move
                    new_board = [list(row) for row in board]
                    new_board[y][x], new_board[ny][nx] = '.', turn
                    new_moves = moves + [f"{turn},{chr(65+x)}{3-y},{chr(65+nx)}{3-ny}"]
                    
                    new_board_tuple = board_to_tuple(new_board)
                    if new_board_tuple not in visited:
                        visited.add(new_board_tuple)
                        next_turn = 'B' if turn == 'w' else 'w'
                        queue.append((new_board, new_moves, next_turn))
    
    return "No"

# Initial board setup
initial_board = [['B', '.', '.', 'B'], ['.', '.', '.', '.'], ['w', 'w', '.', '.']]
result_moves = swap_knights_bfs(initial_board)

print(result_moves)