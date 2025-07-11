import collections

def solve_for_n():
    """
    Finds the value of n for an n x n grid based on a probability problem.
    """
    max_moves = 3
    target_probability = 0.66
    start_pos = (3, 2)  # Corresponds to cell 'c2'

    def is_on_border(pos, n):
        """Checks if a position is on the border of the grid."""
        x, y = pos
        return x == 1 or x == n or y == 1 or y == n

    def get_diagonal_moves(pos, n):
        """Gets all possible cells reachable in one diagonal move."""
        x, y = pos
        moves = set()
        # Up-right
        i, j = x + 1, y + 1
        while i <= n and j <= n:
            moves.add((i, j))
            i += 1
            j += 1
        # Up-left
        i, j = x - 1, y + 1
        while i >= 1 and j <= n:
            moves.add((i, j))
            i -= 1
            j += 1
        # Down-right
        i, j = x + 1, y - 1
        while i <= n and j >= 1:
            moves.add((i, j))
            i += 1
            j -= 1
        # Down-left
        i, j = x - 1, y - 1
        while i >= 1 and j >= 1:
            moves.add((i, j))
            i -= 1
            j -= 1
        return moves

    def get_border_moves(pos, n):
        """Gets all adjacent border cells."""
        x, y = pos
        moves = set()
        if not is_on_border(pos, n):
            return moves
        # Check adjacent cells (up, down, left, right)
        for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nx, ny = x + dx, y + dy
            # The adjacent cell must also be on the grid and on the border
            if 1 <= nx <= n and 1 <= ny <= n and is_on_border((nx, ny), n):
                moves.add((nx, ny))
        return moves

    def count_reachable_cells(n, max_moves, start_pos):
        """
        Counts all cells reachable within max_moves using Breadth-First Search (BFS).
        """
        # A queue for BFS, storing (position, moves_taken)
        q = collections.deque([(start_pos, 0)])
        # A set to store all unique reachable positions
        visited = {start_pos}

        while q:
            pos, moves = q.popleft()

            if moves >= max_moves:
                continue

            # Potential next positions from the current cell
            potential_next_positions = get_diagonal_moves(pos, n)
            if is_on_border(pos, n):
                potential_next_positions.update(get_border_moves(pos, n))

            for next_pos in potential_next_positions:
                if next_pos not in visited:
                    visited.add(next_pos)
                    q.append((next_pos, moves + 1))
        
        return len(visited)

    # We deduced that n must be a multiple of 10.
    # We start checking from n=10 upwards.
    n = 10
    while True:
        total_cells = n * n
        # The target number of reachable cells for the probability to be exactly 0.66
        target_reachable_count = int(target_probability * total_cells)
        
        # Calculate the actual number of reachable cells for the current n
        actual_reachable_count = count_reachable_cells(n, max_moves, start_pos)

        if actual_reachable_count == target_reachable_count:
            print("Solution Found!")
            print(f"For a grid size of n = {n}:")
            print(f"Number of reachable cells within 3 moves: {actual_reachable_count}")
            print(f"Total number of cells in the grid: {n}^2 = {total_cells}")
            print("The probability equation is:")
            print(f"{actual_reachable_count} / {total_cells} = {actual_reachable_count / total_cells}")
            print("\nThis matches the target probability of 66%.")
            print(f"Therefore, the value of n is {n}.")
            return n
        
        # Move to the next multiple of 10
        n += 10
        # Add a reasonable upper limit to prevent an infinite loop
        if n > 200:
            print("No solution found for n up to 200.")
            return None

# Execute the solver
found_n = solve_for_n()
if found_n is not None:
    print(f"<<<{found_n}>>>")
