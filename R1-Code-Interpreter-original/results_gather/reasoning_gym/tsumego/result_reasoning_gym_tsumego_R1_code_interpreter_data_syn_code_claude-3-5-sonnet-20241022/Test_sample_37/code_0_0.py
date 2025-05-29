def print_board(board):
    print("   A B C D E F G H I")
    for i in range(9):
        row = f" {9-i}"
        for j in range(9):
            row += f" {board[i][j]}"
        print(row)

def get_liberties(x, y, board, group=None, visited=None):
    if group is None:
        group = set()
    if visited is None:
        visited = set()
    
    if (x, y) in visited:
        return set()
    
    visited.add((x, y))
    color = board[x][y]
    if color == '.':
        return set()
    
    group.add((x, y))
    liberties = set()
    
    directions = [(-1,0), (1,0), (0,-1), (0,1)]
    for dx, dy in directions:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 9 and 0 <= new_y < 9:
            if board[new_x][new_y] == '.':
                liberties.add((new_x, new_y))
            elif board[new_x][new_y] == color:
                new_liberties = get_liberties(new_x, new_y, board, group, visited)
                liberties.update(new_liberties)
    
    return liberties

# Initialize board
board = [
    ['.', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['X', 'O', 'X', '.', '.', '.', '.', '.', '.'],
    ['O', 'O', '.', 'X', '.', '.', '.', '.', '.'],
    ['X', 'O', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', '.', 'O', '.', '.', '.', '.', '.'],
    ['X', '.', '.', '.', '.', 'O', '.', '.', 'O'],
    ['.', 'X', '.', '.', '.', '.', '.', 'O', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', 'X', '.', '.']
]

# Check liberties of white group at B7
white_liberties = get_liberties(2, 1, board)
print("White group liberties:")
for lib in white_liberties:
    print(f"({chr(lib[1]+65)}{9-lib[0]})")