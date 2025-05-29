def find_captures():
    # Initialize the board (9x9) - corrected position
    board = [
        [".", ".", ".", ".", "O", ".", ".", "O", "."],  # 9
        ["X", ".", ".", ".", "X", "O", "X", ".", "X"],  # 8
        [".", ".", "X", "X", "O", "O", "X", ".", "."],  # 7
        [".", ".", ".", "X", "O", ".", "X", ".", "."],  # 6
        [".", ".", ".", ".", "X", ".", "O", ".", "."],  # 5
        [".", ".", ".", ".", "O", ".", "O", ".", "."],  # 4
        [".", ".", ".", ".", ".", ".", ".", ".", "."],  # 3
        [".", ".", ".", ".", ".", ".", ".", ".", "."],  # 2
        [".", ".", ".", ".", ".", ".", ".", ".", "."]   # 1
    ]

    def is_valid_pos(x, y):
        return 0 <= x < 9 and 0 <= y < 9

    def get_group_liberties(x, y, stone, visited):
        if not is_valid_pos(x, y):
            return set()
        if (x, y) in visited:
            return set()
        if board[x][y] == '.':
            return {(x, y)}
        if board[x][y] != stone:
            return set()

        visited.add((x, y))
        liberties = set()
        for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
            nx, ny = x + dx, y + dy
            liberties.update(get_group_liberties(nx, ny, stone, visited))
        return liberties

    def would_capture(x, y):
        # Place a black stone temporarily
        original = board[x][y]
        board[x][y] = 'X'
        captures = 0
        
        # Check each adjacent position for white groups with no liberties
        for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
            nx, ny = x + dx, y + dy
            if is_valid_pos(nx, ny) and board[nx][ny] == 'O':
                liberties = get_group_liberties(nx, ny, 'O', set())
                if not liberties:
                    # Count stones in captured group
                    visited = set()
                    stack = [(nx, ny)]
                    group_size = 0
                    while stack:
                        cx, cy = stack.pop()
                        if (cx, cy) in visited:
                            continue
                        visited.add((cx, cy))
                        if board[cx][cy] == 'O':
                            group_size += 1
                            for ddx, ddy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                                nnx, nny = cx + ddx, cy + ddy
                                if is_valid_pos(nnx, nny) and board[nnx][nny] == 'O':
                                    stack.append((nnx, nny))
                    captures += group_size

        # Restore the board
        board[x][y] = original
        return captures

    # Check each empty position
    moves = []
    for i in range(9):
        for j in range(9):
            if board[i][j] == '.':
                captures = would_capture(i, j)
                if captures > 0:
                    move = f"{chr(ord('A') + j)}{9-i}"
                    moves.append((move, captures))

    # Sort and print results
    moves.sort(key=lambda x: x[1], reverse=True)
    if moves:
        for move, stones in moves:
            print(f"Move {move}: captures {stones} stones")
    else:
        print("No capturing moves found")

find_captures()