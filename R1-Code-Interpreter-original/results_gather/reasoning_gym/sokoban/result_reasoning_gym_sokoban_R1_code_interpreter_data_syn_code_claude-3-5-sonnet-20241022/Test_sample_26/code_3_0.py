def verify_solution():
    # Initial state
    state = [
        list("+++++++++++"),
        list("++---X---$+"),
        list("++--@-@--++"),
        list("+-X--X-@X+"),
        list("+------$@+"),
        list("+------$*+"),
        list("+++++++++++")
    ]
    
    # Solution to verify
    solution = "LLUURRDLULDRR"
    
    def find_player(board):
        for i in range(len(board)):
            for j in range(len(board[0])):
                if board[i][j] in ['@', '*']:
                    return (i, j)
        return None
    
    def make_move(board, direction):
        board = [row[:] for row in board]
        player = find_player(board)
        if not player:
            return None
            
        px, py = player
        
        # Define move directions
        moves = {
            'U': (-1, 0),
            'D': (1, 0),
            'L': (0, -1),
            'R': (0, 1)
        }
        
        dx, dy = moves[direction]
        new_px, new_py = px + dx, py + dy
        
        # Check if move is valid
        if board[new_px][new_py] == '+':
            return None
            
        # Moving to empty space or goal
        if board[new_px][new_py] in ['-', 'X']:
            # Update player position
            if board[px][py] == '*':
                board[px][py] = 'X'
            else:
                board[px][py] = '-'
                
            if board[new_px][new_py] == 'X':
                board[new_px][new_py] = '*'
            else:
                board[new_px][new_py] = '@'
                
        # Moving a box
        elif board[new_px][new_py] in ['@', '$']:
            box_px, box_py = new_px + dx, new_py + dy
            
            # Check if box can be moved
            if board[box_px][box_py] in ['-', 'X']:
                # Move box
                if board[box_px][box_py] == 'X':
                    board[box_px][box_py] = '$'
                else:
                    board[box_px][box_py] = '@'
                    
                # Move player
                if board[px][py] == '*':
                    board[px][py] = 'X'
                else:
                    board[px][py] = '-'
                    
                if board[new_px][new_py] == '$':
                    board[new_px][new_py] = '*'
                else:
                    board[new_px][new_py] = '@'
            else:
                return None
        else:
            return None
            
        return board
    
    def print_board(board):
        for row in board:
            print(''.join(row))
        print()
    
    # Verify each move
    current_state = state
    print("Initial state:")
    print_board(current_state)
    
    for move in solution:
        new_state = make_move(current_state, move)
        if new_state is None:
            print(f"Invalid move: {move}")
            return False
        current_state = new_state
        print(f"After move {move}:")
        print_board(current_state)
    
    return True

print(verify_solution())