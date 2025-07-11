import sys

# It's recommended to run this with PyPy for better performance, as the search space is large.
# On standard CPython, this script might take a considerable amount of time.

# Set a higher recursion limit for the deep search
sys.setrecursionlimit(2000)

N = 8
grid = [[0] * N for _ in range(N)]
count = 0

def check_row_horizontally(r):
    """Checks a single row for horizontal words with length < 3."""
    c = 0
    while c < N:
        if grid[r][c] == 0:
            start_c = c
            length = 0
            while c < N and grid[r][c] == 0:
                length += 1
                c += 1
            
            is_left_bounded = (start_c == 0) or (grid[r][start_c - 1] == 1)
            is_right_bounded = (c == N) or (grid[r][c] == 1)
            
            if is_left_bounded and is_right_bounded:
                if length < 3:
                    return False
        else:
            c += 1
    return True

def check_final_validity():
    """Performs all final checks on a fully generated grid."""
    
    # 1. Check vertical word lengths
    for c in range(N):
        r = 0
        while r < N:
            if grid[r][c] == 0:
                start_r = r
                length = 0
                while r < N and grid[r][c] == 0:
                    length += 1
                    r += 1
                
                is_top_bounded = (start_r == 0) or (grid[start_r - 1][c] == 1)
                is_bottom_bounded = (r == N) or (grid[r][c] == 1)
                
                if is_top_bounded and is_bottom_bounded:
                    if length < 3:
                        return False
            else:
                r += 1

    # 2. Check connectivity
    white_squares = []
    for r in range(N):
        for c in range(N):
            if grid[r][c] == 0:
                white_squares.append((r, c))

    # Exclude trivial all-black or all-white grids
    if not white_squares or len(white_squares) == N * N:
        return False
    
    q = [white_squares[0]]
    visited = {white_squares[0]}
    head = 0
    while head < len(q):
        r, c = q[head]
        head += 1
        
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nr, nc = r + dr, c + dc
            if 0 <= nr < N and 0 <= nc < N and grid[nr][nc] == 0 and (nr, nc) not in visited:
                visited.add((nr, nc))
                q.append((nr, nc))
                
    if len(visited) != len(white_squares):
        return False

    # 3. Check for cheater squares
    for r in range(N):
        for c in range(N):
            if grid[r][c] == 1:
                def get_color(row, col):
                    if 0 <= row < N and 0 <= col < N:
                        return grid[row][col]
                    return 1 # Treat boundaries as black

                top = get_color(r - 1, c)
                bottom = get_color(r + 1, c)
                left = get_color(r, c - 1)
                right = get_color(r, c + 1)

                # A cheater square is one where flipping it doesn't change the word count.
                # This happens if a new word is created in one dimension while two words are merged in the other.
                if (top == 0 and bottom == 0 and left == 1 and right == 1) or \
                   (left == 0 and right == 0 and top == 1 and bottom == 1):
                    return False
    
    return True

def solve(k):
    """Recursively builds the grid, checking for valid patterns."""
    global count
    
    # Base case: if the top half of the grid is filled
    if k == N * N // 2:
        if check_final_validity():
            count += 1
        return

    r = k // N
    c = k % N
    sr, sc = N - 1 - r, N - 1 - c

    # Branch 1: Try placing a white square
    grid[r][c] = 0
    grid[sr][sc] = 0
    
    # Prune if a row is complete and invalid
    should_prune = (c == N - 1) and not (check_row_horizontally(r) and check_row_horizontally(sr))
    
    if not should_prune:
        solve(k + 1)

    # Branch 2: Try placing a black square
    grid[r][c] = 1
    grid[sr][sc] = 1

    should_prune = (c == N - 1) and not (check_row_horizontally(r) and check_row_horizontally(sr))

    if not should_prune:
        solve(k + 1)


if __name__ == '__main__':
    solve(0)
    print(f"The number of possible 8x8 crossword grids is:")
    print(count)