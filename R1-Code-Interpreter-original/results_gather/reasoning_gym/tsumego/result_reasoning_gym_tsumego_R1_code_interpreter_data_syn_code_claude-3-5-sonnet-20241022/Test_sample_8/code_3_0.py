def analyze_board():
    # Create the board as a 2D list (12x12)
    board = [
        ['.', '.', 'O', '.', '.', '.', '.', 'X', '.', '.', '.', '.'],  # 12
        ['.', '.', '.', '.', '.', '.', 'X', 'O', 'X', '.', '.', '.'],  # 11
        ['.', '.', '.', '.', '.', 'X', 'O', 'O', '.', '.', '.', '.'],  # 10
        ['.', '.', '.', '.', '.', '.', 'X', 'O', 'X', 'O', '.', '.'],  # 9
        ['.', '.', '.', '.', '.', '.', '.', 'X', '.', '.', 'O', '.'],  # 8
        ['.', 'X', 'O', '.', 'O', '.', '.', '.', '.', '.', '.', '.'],  # 7
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 6
        ['O', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 5
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 4
        ['.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.', '.'],  # 3
        ['.', '.', '.', 'X', '.', '.', '.', '.', 'O', '.', '.', '.'],  # 2
        ['.', '.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.']   # 1
    ]

    def is_valid(x, y):
        return 0 <= x < 12 and 0 <= y < 12

    def get_group(x, y, color, visited=None):
        if visited is None:
            visited = set()
        
        if not is_valid(x, y) or board[x][y] != color or (x, y) in visited:
            return set()
        
        visited.add((x, y))
        group = {(x, y)}
        
        for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            new_x, new_y = x + dx, y + dy
            group.update(get_group(new_x, new_y, color, visited))
        
        return group

    def count_liberties(group):
        liberties = set()
        for x, y in group:
            for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                new_x, new_y = x + dx, y + dy
                if is_valid(new_x, new_y) and board[new_x][new_y] == '.':
                    liberties.add((new_x, new_y))
        return liberties

    # Find all white groups and their liberties
    white_groups = []
    visited = set()
    
    for i in range(12):
        for j in range(12):
            if board[i][j] == 'O' and (i, j) not in visited:
                group = get_group(i, j, 'O')
                visited.update(group)
                liberties = count_liberties(group)
                white_groups.append((group, liberties))

    # Check for potential captures
    capturing_moves = []
    for group, liberties in white_groups:
        if len(liberties) == 1:
            liberty = list(liberties)[0]
            # Convert to board coordinates (A-L, 12-1)
            move = f"{chr(ord('A') + liberty[1])}{12-liberty[0]}"
            capturing_moves.append((move, len(group)))

    print(f"White groups found: {len(white_groups)}")
    for i, (group, liberties) in enumerate(white_groups):
        print(f"Group {i+1}: {len(group)} stones, {len(liberties)} liberties")
    print(f"Potential capturing moves: {capturing_moves}")

analyze_board()