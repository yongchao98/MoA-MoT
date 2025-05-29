from collections import deque

def get_valid_moves(board, pos):
    x, y = ord(pos[0]) - ord('A'), int(pos[1]) - 1
    moves = []
    for dx, dy in [(2,1), (2,-1), (-2,1), (-2,-1), (1,2), (1,-2), (-1,2), (-1,-2)]:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 4 and 0 <= new_y < 3:  # 4x3 board
            new_pos = f"{chr(ord('A') + new_x)}{new_y + 1}"
            if board[new_pos] == '.':
                moves.append(new_pos)
    return moves

def board_to_string(board):
    return ''.join(board[pos] for pos in sorted(board.keys()))

def solve_knight_swap():
    # Initialize board
    board = {}
    for x in 'ABCD':
        for y in range(1, 4):
            board[f"{x}{y}"] = '.'
    
    # Set initial positions
    board['B3'] = 'B'
    board['C2'] = 'w'
    board['D2'] = 'w'
    board['C1'] = 'B'
    
    # Target positions for white knights
    target_positions = {'B3', 'C1'}
    
    start_state = (board, True)  # (board, is_white_turn)
    queue = deque([(start_state, [])])  # (state, moves)
    visited = set()
    
    while queue:
        (current_board, is_white_turn), moves = queue.popleft()
        
        # Create unique state key
        state_key = (board_to_string(current_board), is_white_turn)
        if state_key in visited:
            continue
        visited.add(state_key)
        
        # Check if target reached
        white_positions = {pos for pos, piece in current_board.items() if piece == 'w'}
        if white_positions == target_positions and len(moves) > 0:
            return moves
        
        # Generate moves for current player
        piece = 'w' if is_white_turn else 'B'
        positions = [pos for pos, p in current_board.items() if p == piece]
        
        for pos in positions:
            for new_pos in get_valid_moves(current_board, pos):
                # Make move
                new_board = current_board.copy()
                new_board[pos] = '.'
                new_board[new_pos] = piece
                
                # Add move
                move = f"{piece.lower()},{pos},{new_pos}"
                queue.append(((new_board, not is_white_turn), moves + [move]))
    
    return None

# Run solution
result = solve_knight_swap()
if result:
    print(result)
else:
    print("No")