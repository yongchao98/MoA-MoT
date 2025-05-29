def create_board():
    board = [
        ['X', '.', '.', '.', 'X', '.', '.', '.', '.'],  # 9
        ['.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 8
        ['.', 'O', '.', 'O', '.', '.', 'X', '.', '.'],  # 7
        ['.', '.', '.', 'X', '.', '.', '.', '.', 'O'],  # 6
        ['O', '.', 'X', 'O', 'X', '.', '.', '.', '.'],  # 5
        ['.', 'X', 'O', 'O', '.', 'O', '.', '.', '.'],  # 4
        ['.', '.', 'X', 'O', 'X', '.', '.', '.', '.'],  # 3
        ['.', '.', '.', 'X', '.', '.', '.', '.', '.'],  # 2
        ['.', 'O', '.', 'O', '.', '.', 'X', '.', '.']   # 1
    ]
    return board

def get_liberties(board, x, y, visited):
    if x < 0 or x >= 9 or y < 0 or y >= 9:
        return 0
    if board[y][x] == '.':
        return 1
    if (x, y) in visited:
        return 0
    if board[y][x] == 'O':
        visited.add((x, y))
        liberties = 0
        # Check all four directions
        directions = [(0, 1), (0, -1), (1, 0), (-1, 0)]
        for dx, dy in directions:
            liberties += get_liberties(board, x + dx, y + dy, visited)
        return liberties
    return 0

def find_capturing_moves():
    board = create_board()
    potential_moves = []
    
    # Check each empty position
    for y in range(9):
        for x in range(9):
            if board[y][x] == '.':
                # Try each empty position
                board[y][x] = 'X'  # Temporarily place black stone
                
                # Check if this move captures any white groups
                for ny in range(9):
                    for nx in range(9):
                        if board[ny][nx] == 'O':
                            liberties = get_liberties(board, nx, ny, set())
                            if liberties == 0:
                                move = chr(ord('A') + x) + str(9 - y)
                                potential_moves.append((move, liberties))
                
                board[y][x] = '.'  # Reset the position
    
    print("Potential capturing moves:", potential_moves)

find_capturing_moves()