def solve_rush_hour():
    # Initial state analysis
    # AA is in row 2 (0-based), needs to reach right edge
    # Blocking cars: G, I, J, K
    # Need to clear the path in row 2
    
    def is_path_clear(moves):
        board = [
            list('BB.HJK'),
            list('..xHJK'),
            list('AAGIJK'),
            list('.FGI..'),
            list('.FCC..'),
            list('.DDEE.')
        ]
        
        # Apply moves
        for move in moves:
            car = move[0]
            direction = 1 if '+' in move else -1
            spaces = int(move[-1])
            
            # Find car positions
            positions = []
            for i in range(6):
                for j in range(6):
                    if board[i][j] == car:
                        positions.append((i, j))
            
            # Move car
            is_horizontal = positions[0][0] == positions[1][0]
            # Clear current positions
            for pos in positions:
                board[pos[0]][pos[1]] = '.'
            
            # Set new positions
            if is_horizontal:
                for pos in positions:
                    board[pos[0]][pos[1] + direction * spaces] = car
            else:
                for pos in positions:
                    board[pos[0] + direction * spaces][pos[1]] = car
        
        # Check if path is clear for AA to exit
        row = board[2]
        start_pos = None
        for j in range(6):
            if row[j] == 'A':
                start_pos = j
                break
        
        if start_pos is None:
            return False
            
        # Check if path to exit is clear
        for j in range(start_pos + 2, 6):
            if row[j] != '.':
                return False
        return True

    # Try the most direct sequence to clear the path
    moves = ['J-1', 'K-1', 'I-1', 'G-1', 'A+3']
    
    if is_path_clear(moves):
        return moves
    
    return None

solution = solve_rush_hour()
print('<<<' + ' '.join(solution) + '>>>')