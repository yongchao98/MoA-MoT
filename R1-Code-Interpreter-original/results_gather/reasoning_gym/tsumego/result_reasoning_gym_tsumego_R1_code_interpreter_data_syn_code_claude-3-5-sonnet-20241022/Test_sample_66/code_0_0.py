def print_board(board):
    print("   A B C D E F G H I J K")
    for i in range(11):
        row = 11 - i
        print(f"{row:2d} {' '.join(board[i])}")

def find_liberties(board, x, y, color, visited):
    if x < 0 or x >= 11 or y < 0 or y >= 11:
        return set()
    if board[y][x] == '.':
        return {(x, y)}
    if board[y][x] != color or (x, y) in visited:
        return set()
    
    visited.add((x, y))
    liberties = set()
    for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
        liberties.update(find_liberties(board, x + dx, y + dy, color, visited))
    return liberties

def analyze_position():
    board = [
        ['.', '.', '.', '.', '.', '.', 'O', 'X', '.', '.', '.'],  # 11
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 10
        ['.', '.', '.', 'O', '.', '.', '.', '.', '.', '.', '.'],  # 9
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 8
        ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'O', '.'],  # 7
        ['O', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],  # 6
        ['.', '.', '.', '.', '.', '.', '.', '.', 'X', 'O', 'X'],  # 5
        ['.', '.', '.', 'X', '.', '.', '.', 'X', 'O', 'O', '.'],  # 4
        ['.', '.', 'X', '.', '.', '.', '.', '.', 'X', 'O', 'X'],  # 3
        ['X', '.', '.', '.', '.', '.', '.', '.', 'O', 'X', '.'],  # 2
        ['.', '.', '.', '.', '.', 'X', '.', '.', 'X', '.', '.'],  # 1
    ]

    # Find all white groups and their liberties
    white_groups = []
    visited = set()
    
    for y in range(11):
        for x in range(11):
            if board[y][x] == 'O' and (x, y) not in visited:
                liberties = find_liberties(board, x, y, 'O', visited)
                if liberties:
                    white_groups.append({
                        'position': (x, y),
                        'liberties': liberties
                    })

    # Print analysis
    print("White groups and their liberties:")
    cols = 'ABCDEFGHIJK'
    for group in white_groups:
        x, y = group['position']
        pos = f"{cols[x]}{11-y}"
        lib_positions = [f"{cols[x]}{11-y}" for x, y in group['liberties']]
        print(f"Group at {pos} has liberties at: {', '.join(lib_positions)}")

analyze_position()