from collections import deque

def get_knight_moves(pos):
    row, col = pos
    moves = [
        (row+2, col+1), (row+2, col-1),
        (row-2, col+1), (row-2, col-1),
        (row+1, col+2), (row+1, col-2),
        (row-1, col+2), (row-1, col-2)
    ]
    return [(r,c) for r,c in moves if 0 <= r < 3 and 0 <= c < 3]

def pos_to_notation(pos):
    return f"{chr(65+pos[1])}{3-pos[0]}"

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def find_pieces(board, piece):
    return [(i,j) for i in range(3) for j in range(3) if board[i][j] == piece]

def is_target_reached(board):
    white_positions = set((r,c) for r in range(3) for c in range(3) if board[r][c] == 'w')
    black_positions = set((r,c) for r in range(3) for c in range(3) if board[r][c] == 'B')
    return white_positions == {(2,0), (2,2)} and black_positions == {(0,2), (1,0)}

def solve_with_depth_limit(max_depth=12):  # Set reasonable max depth
    initial_board = [
        ['.', '.', 'w'],  # row 3
        ['w', '.', '.'],  # row 2
        ['B', '.', 'B']   # row 1
    ]
    
    def board_to_string(board):
        return ''.join(''.join(row) for row in board)
    
    def make_move(board, start, end):
        new_board = [row[:] for row in board]
        piece = new_board[start[0]][start[1]]
        new_board[start[0]][start[1]] = '.'
        new_board[end[0]][end[1]] = piece
        return new_board
    
    def dfs(board, moves, depth, is_white_turn, seen):
        if depth >= max_depth:
            return None
            
        if is_target_reached(board):
            return moves
            
        current_piece = 'w' if is_white_turn else 'B'
        pieces = find_pieces(board, current_piece)
        
        # Sort moves by proximity to target positions
        possible_moves = []
        for start_pos in pieces:
            for end_pos in get_knight_moves(start_pos):
                if board[end_pos[0]][end_pos[1]] == '.':
                    # Calculate score based on target positions
                    score = 0
                    if current_piece == 'w':
                        score = -manhattan_distance(end_pos, (2,0)) - manhattan_distance(end_pos, (2,2))
                    else:
                        score = -manhattan_distance(end_pos, (0,2)) - manhattan_distance(end_pos, (1,0))
                    possible_moves.append((score, start_pos, end_pos))
        
        # Try moves in order of best score
        for _, start_pos, end_pos in sorted(possible_moves, reverse=True):
            new_board = make_move(board, start_pos, end_pos)
            board_str = board_to_string(new_board)
            
            if board_str not in seen:
                seen.add(board_str)
                move = f"{current_piece},{pos_to_notation(start_pos)},{pos_to_notation(end_pos)}"
                result = dfs(new_board, moves + [move], depth + 1, not is_white_turn, seen)
                if result:
                    return result
                seen.remove(board_str)
                
        return None

    seen = {board_to_string(initial_board)}
    return dfs(initial_board, [], 0, True, seen)

solution = solve_with_depth_limit()
if solution is None:
    print("<<<No>>>")
else:
    print(f"<<<{solution}>>>")