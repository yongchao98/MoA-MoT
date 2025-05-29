def find_captures():
    # Board representation (0 = empty, 1 = black, 2 = white)
    board = [
        [0, 0, 0, 0, 2, 0, 0, 2, 0],  # 9
        [1, 0, 0, 0, 1, 2, 1, 0, 1],  # 8
        [0, 0, 1, 1, 2, 2, 1, 0, 0],  # 7
        [0, 0, 0, 1, 2, 0, 1, 0, 0],  # 6
        [0, 0, 0, 0, 1, 0, 2, 0, 0],  # 5
        [0, 0, 0, 0, 2, 0, 2, 0, 0],  # 4
        [0, 0, 0, 0, 0, 0, 0, 0, 0],  # 3
        [0, 0, 0, 0, 0, 0, 0, 0, 0],  # 2
        [0, 0, 0, 0, 0, 0, 0, 0, 0]   # 1
    ]

    def has_liberty(x, y, visited=None):
        if visited is None:
            visited = set()
            
        if x < 0 or x >= 9 or y < 0 or y >= 9:
            return False
            
        if board[x][y] == 0:  # Empty spot = liberty
            return True
            
        if (x, y) in visited:
            return False
            
        if board[x][y] != 2:  # Not white stone
            return False
            
        visited.add((x, y))
        
        # Check all adjacent positions
        for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
            new_x, new_y = x + dx, y + dy
            if has_liberty(new_x, new_y, visited):
                return True
                
        return False

    def count_captured_stones(x, y, visited=None):
        if visited is None:
            visited = set()
            
        if x < 0 or x >= 9 or y < 0 or y >= 9:
            return 0
            
        if (x, y) in visited or board[x][y] != 2:
            return 0
            
        visited.add((x, y))
        count = 1
        
        for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
            new_x, new_y = x + dx, y + dy
            count += count_captured_stones(new_x, new_y, visited)
            
        return count

    # Try each empty position
    captures = []
    for i in range(9):
        for j in range(9):
            if board[i][j] == 0:  # Empty position
                # Temporarily place black stone
                board[i][j] = 1
                
                # Check adjacent positions for white groups without liberties
                total_captured = 0
                checked = set()
                for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                    new_x, new_y = i + dx, j + dy
                    if (0 <= new_x < 9 and 0 <= new_y < 9 and 
                        board[new_x][new_y] == 2 and 
                        (new_x, new_y) not in checked):
                        if not has_liberty(new_x, new_y):
                            total_captured += count_captured_stones(new_x, new_y)
                            
                # Remove temporary stone
                board[i][j] = 0
                
                if total_captured > 0:
                    move = f"{chr(ord('A') + j)}{9-i}"
                    captures.append((move, total_captured))
                    print(f"Move {move}: captures {total_captured} stones")

find_captures()