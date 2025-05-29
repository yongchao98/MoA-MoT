from collections import deque

def get_knight_moves(x, y):
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(nx,ny) for nx,ny in moves if 0 <= nx < 4 and 0 <= ny < 3]

def encode_state(board, white_turn):
    return (''.join(''.join(row) for row in board), white_turn)

def is_target_state(board):
    # Target: White at B3,C3 (where black started), Black at A2,D2 (where white started)
    return (board[0][1] == 'w' and board[0][2] == 'w' and  # B3,C3
            board[1][0] == 'B' and board[1][3] == 'B')     # A2,D2

def solve():
    # Initial board: Black at B3,C3; White at A2,D2
    initial_board = [
        ['.', 'B', 'B', '.'],  # row 3 (index 0)
        ['w', '.', '.', 'w'],  # row 2 (index 1)
        ['.', '.', '.', '.']   # row 1 (index 2)
    ]
    
    queue = deque([(initial_board, [], True)])  # board, moves, white_turn
    visited = {encode_state(initial_board, True)}
    
    while queue:
        board, moves, white_turn = queue.popleft()
        piece = 'w' if white_turn else 'B'
        
        # Find current player's knights
        for y in range(3):
            for x in range(4):
                if board[y][x] == piece:
                    # Try all possible moves for this knight
                    for new_x, new_y in get_knight_moves(x, y):
                        if board[new_y][new_x] == '.':  # Empty square
                            # Make move
                            new_board = [row[:] for row in board]
                            new_board[y][x] = '.'
                            new_board[new_y][new_x] = piece
                            
                            # Generate move notation
                            from_pos = f"{chr(x + ord('A'))}{3-y}"
                            to_pos = f"{chr(new_x + ord('A'))}{3-new_y}"
                            new_moves = moves + [f"{piece},{from_pos},{to_pos}"]
                            
                            # Check if this is a new state
                            state = encode_state(new_board, not white_turn)
                            if state not in visited:
                                visited.add(state)
                                
                                # Check if we reached target state
                                if is_target_state(new_board):
                                    return new_moves
                                
                                queue.append((new_board, new_moves, not white_turn))
    return "No"

solution = solve()
print(solution)