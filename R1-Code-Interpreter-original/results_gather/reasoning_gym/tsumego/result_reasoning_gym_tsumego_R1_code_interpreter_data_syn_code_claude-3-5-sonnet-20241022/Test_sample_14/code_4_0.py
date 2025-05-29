def analyze_liberties(move_x, move_y):
    # Initialize the board
    board = [
        ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', 'O', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.', '.', 'O'],
        ['.', '.', '.', '.', '.', '.', '.', '.', 'X'],
        ['.', '.', '.', 'X', '.', '.', '.', '.', '.'],
        ['O', '.', 'X', 'O', 'X', '.', '.', '.', 'X'],
        ['O', 'X', 'O', 'O', '.', '.', '.', '.', '.'],
        ['.', 'X', 'X', 'O', 'X', 'X', '.', '.', '.'],
        ['X', 'X', '.', 'X', '.', 'O', '.', '.', '.']
    ]
    
    # Make the move
    board[move_y][move_x] = 'X'
    
    def count_liberties(x, y, color, visited=None):
        if visited is None:
            visited = set()
            
        if x < 0 or x >= 9 or y < 0 or y >= 9:
            return set()
            
        if board[y][x] != color:
            if board[y][x] == '.':
                return {(x, y)}
            return set()
            
        if (x, y) in visited:
            return set()
            
        visited.add((x, y))
        liberties = set()
        
        for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            new_x, new_y = x + dx, y + dy
            if 0 <= new_x < 9 and 0 <= new_y < 9:
                if board[new_y][new_x] == '.':
                    liberties.add((new_x, new_y))
                elif board[new_y][new_x] == color:
                    liberties.update(count_liberties(new_x, new_y, color, visited))
        
        return liberties

    # Check white group starting from A3 (0,6)
    liberties = count_liberties(0, 6, 'O')
    print(f"After move at {chr(move_x+65)}{9-move_y}, white group liberties:", 
          [(chr(x+65), 9-y) for x, y in liberties])
    return len(liberties)

# Test the three critical points
print("Testing A5:")
liberties_a5 = analyze_liberties(0, 4)
print("\nTesting B6:")
liberties_b6 = analyze_liberties(1, 3)
print("\nTesting A8:")
liberties_a8 = analyze_liberties(0, 1)