def analyze_capture_sequence(board):
    def count_liberties(x, y, color, visited=None):
        if visited is None:
            visited = set()
        if (x, y) in visited:
            return set()
        if x < 0 or x >= 9 or y < 0 or y >= 9:
            return set()
        if board[y][x] != color:
            return {(x,y)} if board[y][x] == '.' else set()
            
        visited.add((x, y))
        liberties = set()
        for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
            liberties.update(count_liberties(x + dx, y + dy, color, visited))
        return liberties

    # Make move at G5
    test_board = [row[:] for row in board]
    test_board[4][6] = 'X'  # G5
    
    # Check upper white group liberties
    upper_group_stones = [(6,2), (7,2), (7,1), (7,3)]  # G7, H7, H8, H6
    upper_liberties = set()
    visited = set()
    
    for x, y in upper_group_stones:
        upper_liberties.update(count_liberties(x, y, 'O', visited))
    
    print("After Black G5:")
    print("Upper white group liberties:", len(upper_liberties))
    print("Liberty points:", [(chr(65+x), 9-y) for x,y in upper_liberties])
    
    # Check if white can escape
    can_escape = False
    for x, y in upper_liberties:
        # Check if white can make new liberties by playing at liberty point
        test_board[y][x] = 'O'
        new_liberties = count_liberties(x, y, 'O')
        if len(new_liberties) > 1:
            can_escape = True
        test_board[y][x] = '.'
    
    print("Can white escape?", can_escape)

# Initialize board
board = [
    ['.', '.', 'O', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', 'X', '.', '.', '.', 'X', 'O', 'X'],
    ['X', 'X', '.', '.', '.', 'X', 'O', 'O', '.'],
    ['.', '.', '.', 'X', 'X', '.', 'X', 'O', 'X'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', 'O', '.', '.', 'O', '.', '.', '.'],
    ['X', '.', '.', '.', '.', 'O', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.']
]

analyze_capture_sequence(board)