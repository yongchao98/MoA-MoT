from collections import deque

def get_knight_moves(pos):
    row, col = pos
    moves = [
        (row+2, col+1), (row+2, col-1),
        (row-2, col+1), (row-2, col-1),
        (row+1, col+2), (row+1, col-2),
        (row-1, col+2), (row-1, col-2)
    ]
    return [(r, c) for r, c in moves if 0 <= r < 3 and 0 <= c < 4]

def pos_to_notation(pos):
    return f"{chr(65+pos[1])}{3-pos[0]}"

def get_pieces_positions(board):
    whites = []
    blacks = []
    empty = []
    for i in range(3):
        for j in range(4):
            if board[i][j] == 'w':
                whites.append((i, j))
            elif board[i][j] == 'B':
                blacks.append((i, j))
            elif board[i][j] == '.':
                empty.append((i, j))
    return whites, blacks, empty

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_target_reached(whites, blacks):
    # Target positions for white knights: C3 and A2 (original black positions)
    white_target = {(0, 2), (1, 0)}
    # Target positions for black knights: B3 and C1 (original white positions)
    black_target = {(0, 1), (2, 2)}
    
    return (set((r, c) for r, c in whites) == white_target and 
            set((r, c) for r, c in blacks) == black_target)

def find_shortest_solution():
    initial_board = [
        ['.', 'w', 'B', ' '],  # row 3
        ['B', ' ', ' ', '.'],  # row 2
        [' ', '.', 'w', ' ']   # row 1
    ]
    
    queue = deque([(initial_board, [], True)])
    visited = {board_to_string(initial_board)}
    max_depth = 10  # Limit search depth
    
    while queue:
        current_board, moves, is_white_turn = queue.popleft()
        
        if len(moves) >= max_depth:
            continue
            
        whites, blacks, empty = get_pieces_positions(current_board)
        
        if is_target_reached(whites, blacks):
            return moves
        
        pieces = whites if is_white_turn else blacks
        piece_type = 'w' if is_white_turn else 'B'
        
        for piece in pieces:
            for new_pos in get_knight_moves(piece):
                if new_pos in [(r, c) for r, c in empty]:
                    new_board = [row[:] for row in current_board]
                    new_board[piece[0]][piece[1]] = '.'
                    new_board[new_pos[0]][new_pos[1]] = piece_type
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        new_moves = moves + [f"{piece_type},{pos_to_notation(piece)},{pos_to_notation(new_pos)}"]
                        queue.append((new_board, new_moves, not is_white_turn))
    
    return None

solution = find_shortest_solution()
if solution:
    print(solution)
else:
    print("No")