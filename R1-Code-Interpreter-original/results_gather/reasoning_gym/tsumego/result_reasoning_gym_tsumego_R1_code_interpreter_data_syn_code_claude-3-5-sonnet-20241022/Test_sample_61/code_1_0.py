def check_capture_after_move():
    board = [['.'] * 12 for _ in range(12)]
    # Place black stones (X)
    black_stones = [(11,1), (10,6), (10,9), (8,9), (7,7), (7,9), (6,6), (6,10), 
                    (5,7), (5,9), (4,8), (3,11), (2,7), (5,6)]  # Added G5
    # Place white stones (O)
    white_stones = [(9,7), (7,6), (6,7), (6,8), (6,9), (5,8), (4,1), (2,1), (2,3)]
    
    for x, y in black_stones:
        board[x][y] = 'X'
    for x, y in white_stones:
        board[x][y] = 'O'

    def count_liberties(x, y, color, visited=None):
        if visited is None:
            visited = set()
        if (x, y) in visited or x < 0 or x >= 12 or y < 0 or y >= 12:
            return 0
        if board[x][y] != color:
            return 1 if board[x][y] == '.' else 0
        
        visited.add((x,y))
        liberties = 0
        for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
            liberties += count_liberties(x+dx, y+dy, color, visited)
        return liberties

    # Check liberties of white group after G5
    white_group_liberties = count_liberties(6, 7, 'O')
    print(f"White group liberties after Black plays G5: {white_group_liberties}")
    
    # Check if white can make eyes
    potential_eyes = [(7,6), (6,8)]
    eye_spaces = 0
    for x, y in potential_eyes:
        surrounded = True
        for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
            nx, ny = x + dx, y + dy
            if 0 <= nx < 12 and 0 <= ny < 12:
                if board[nx][ny] != 'O' and board[nx][ny] != 'X':
                    surrounded = False
        if surrounded:
            eye_spaces += 1
    print(f"Potential eye spaces: {eye_spaces}")

check_capture_after_move()