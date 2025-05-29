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
    return [(r,c) for r,c in moves if 0 <= r < 3 and 0 <= c < 4]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_target_reached(board):
    # Check if white knights are at original black positions and vice versa
    target_whites = {(1,3), (2,2)}  # Original black positions
    target_blacks = {(2,1), (0,1)}  # Original white positions
    
    current_whites = set()
    current_blacks = set()
    for i in range(3):
        for j in range(4):
            if board[i][j] == 'w':
                current_whites.add((i,j))
            elif board[i][j] == 'B':
                current_blacks.add((i,j))
    
    return current_whites == target_whites and current_blacks == target_blacks

def solve_knight_swap():
    # Initial board
    initial_board = [
        ['.', 'w', '.', '.'],
        ['.', '.', '.', 'B'],
        ['w', '.', 'B', '.']
    ]
    
    # BFS queue: (board, moves_list, is_black_turn)
    queue = deque([(initial_board, [], True)])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves, is_black_turn = queue.popleft()
        
        if is_target_reached(current_board):
            return moves
        
        # Find current player's knights
        piece = 'B' if is_black_turn else 'w'
        knights = []
        for i in range(3):
            for j in range(4):
                if current_board[i][j] == piece:
                    knights.append((i,j))
        
        # Try all possible moves for each knight
        for knight in knights:
            for new_pos in get_knight_moves(knight):
                if current_board[new_pos[0]][new_pos[1]] == '.':
                    # Make move
                    new_board = copy.deepcopy(current_board)
                    new_board[knight[0]][knight[1]] = '.'
                    new_board[new_pos[0]][new_pos[1]] = piece
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        # Convert coordinates to chess notation
                        from_sq = f"{chr(knight[1]+65)}{3-knight[0]}"
                        to_sq = f"{chr(new_pos[1]+65)}{3-new_pos[0]}"
                        new_moves = moves + [f"{piece},{from_sq},{to_sq}"]
                        queue.append((new_board, new_moves, not is_black_turn))
    
    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")