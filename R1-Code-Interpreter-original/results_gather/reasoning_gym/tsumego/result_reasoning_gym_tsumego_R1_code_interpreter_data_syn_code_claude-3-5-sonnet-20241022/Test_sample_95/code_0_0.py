def print_board(board):
    print("   A B C D E F G H I")
    for i in range(9):
        row = f" {9-i}"
        for j in range(9):
            row += f" {board[i][j]}"
        print(row)

def get_liberties(board, x, y, color, visited=None):
    if visited is None:
        visited = set()
    
    if x < 0 or x >= 9 or y < 0 or y >= 9:
        return set()
    
    if board[x][y] == '.':
        return {(x, y)}
    
    if board[x][y] != color or (x, y) in visited:
        return set()
        
    visited.add((x, y))
    liberties = set()
    
    # Check all adjacent positions
    for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
        new_x, new_y = x + dx, y + dy
        liberties.update(get_liberties(board, new_x, new_y, color, visited))
    
    return liberties

def find_captures(board):
    potential_moves = []
    
    # Convert string representation to 2D array
    board = [list(row) for row in [
        ".........",
        "...O...O.",
        "X...XOX.X",
        "..XXOOX..",
        "...XO.X..",
        "....X.O..",
        "....O.O..",
        ".........",
        "........."
    ]]
    
    # Check each empty position
    for i in range(9):
        for j in range(9):
            if board[i][j] == '.':
                # Try placing a black stone here
                board[i][j] = 'X'
                
                # Check if this move captures any white groups
                captured = 0
                for x in range(9):
                    for y in range(9):
                        if board[x][y] == 'O':
                            liberties = get_liberties(board, x, y, 'O')
                            if not liberties:
                                # Count stones in this group
                                group_size = 1
                                temp_board = [row[:] for row in board]
                                stack = [(x, y)]
                                temp_board[x][y] = '.'
                                while stack:
                                    cx, cy = stack.pop()
                                    for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                                        nx, ny = cx + dx, cy + dy
                                        if 0 <= nx < 9 and 0 <= ny < 9 and temp_board[nx][ny] == 'O':
                                            group_size += 1
                                            temp_board[nx][ny] = '.'
                                            stack.append((nx, ny))
                                captured += group_size
                
                if captured > 0:
                    potential_moves.append((chr(ord('A') + j) + str(9-i), captured))
                
                # Reset the position
                board[i][j] = '.'
    
    # Sort moves by number of captures
    potential_moves.sort(key=lambda x: x[1], reverse=True)
    
    print("\nPotential capturing moves:")
    for move, stones in potential_moves:
        print(f"Move {move}: captures {stones} stones")

find_captures(None)