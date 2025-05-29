from collections import deque
import copy

def print_board(board):
    for i, row in enumerate(board):
        print(f"{3-i} |", end=" ")
        print(" | ".join(cell if cell != '.' else ' ' for cell in row), end=" ")
        print("|")
        print("  " + "-" * 19)
    print("    A   B   C   D")

def get_knight_moves(pos):
    row, col = pos
    moves = [
        (row+2, col+1), (row+2, col-1),
        (row-2, col+1), (row-2, col-1),
        (row+1, col+2), (row+1, col-2),
        (row-1, col+2), (row-1, col-2)
    ]
    return [(r,c) for r,c in moves if 0 <= r < 4 and 0 <= c < 4]

def is_target_reached(board):
    # Target positions:
    # White knights should be at original black positions (B2, C3)
    # Black knights should be at original white positions (B1, D1)
    white_positions = set((r,c) for r in range(4) for c in range(4) if board[r][c] == 'w')
    black_positions = set((r,c) for r in range(4) for c in range(4) if board[r][c] == 'B')
    
    target_white = {(1,1), (0,2)}  # B2, C3
    target_black = {(2,1), (2,3)}  # B1, D1
    
    return white_positions == target_white and black_positions == target_black

def solve_knight_swap():
    initial_board = [
        ['.', '.', 'B', '.'],
        ['.', 'B', '.', '.'],
        ['.', 'w', '.', 'w'],
        ['.', '.', '.', '.']
    ]
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    seen = {str(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        if is_target_reached(board):
            return moves
            
        piece = 'w' if is_white_turn else 'B'
        positions = [(r,c) for r in range(4) for c in range(4) if board[r][c] == piece]
        
        for pos in positions:
            for new_pos in get_knight_moves(pos):
                if board[new_pos[0]][new_pos[1]] == '.':
                    new_board = copy.deepcopy(board)
                    new_board[pos[0]][pos[1]] = '.'
                    new_board[new_pos[0]][new_pos[1]] = piece
                    
                    board_str = str(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        from_sq = f"{chr(pos[1]+65)}{3-pos[0]+1}"
                        to_sq = f"{chr(new_pos[1]+65)}{3-new_pos[0]+1}"
                        move = f"{piece},{from_sq},{to_sq}"
                        queue.append((new_board, moves + [move], not is_white_turn))
    
    return "No"

# Test the solution
solution = solve_knight_swap()
print(solution)