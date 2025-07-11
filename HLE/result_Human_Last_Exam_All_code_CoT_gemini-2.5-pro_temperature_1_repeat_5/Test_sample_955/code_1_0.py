import collections

def find_grid_size():
    """
    This function simulates the object's movement on an n x n grid to find
    the value of n that satisfies the given probability condition.

    The probability of reaching a cell within 3 moves is 66%, so:
    Reachable Cells / (n * n) = 0.66 = 33 / 50.
    This implies n*n must be a multiple of 50, so n must be a multiple of 10.
    We will test n = 10, 20, 30, and so on.
    """
    n = 10
    while True:
        # Starting position 'c2' corresponds to coordinates (3, 2)
        start_pos = (3, 2)
        
        # The queue for the BFS algorithm will store ((x, y), distance)
        queue = collections.deque([(start_pos, 0)])
        visited = {start_pos}

        # Perform the BFS to find all cells reachable within 3 moves
        while queue:
            (x, y), dist = queue.popleft()

            # Stop exploring from a cell if it's already 3 moves away
            if dist >= 3:
                continue

            # 1. Generate diagonal moves
            # The signs determine the direction: UR, UL, DR, DL
            for dx_sign in [-1, 1]:
                for dy_sign in [-1, 1]:
                    # Explore along the entire diagonal line
                    for i in range(1, n + 1):
                        nx, ny = x + i * dx_sign, y + i * dy_sign
                        if 1 <= nx <= n and 1 <= ny <= n:
                            if (nx, ny) not in visited:
                                visited.add((nx, ny))
                                queue.append(((nx, ny), dist + 1))
                        else:
                            # Stop when the diagonal goes off the grid
                            break
            
            # 2. Generate border moves
            is_on_border = (x == 1 or x == n or y == 1 or y == n)
            if is_on_border:
                # Explore adjacent cells (N, S, E, W)
                for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                    nx, ny = x + dx, y + dy
                    if 1 <= nx <= n and 1 <= ny <= n:
                        # A border move requires the destination to also be on a border
                        is_neighbor_on_border = (nx == 1 or nx == n or ny == 1 or ny == n)
                        if is_neighbor_on_border:
                            if (nx, ny) not in visited:
                                visited.add((nx, ny))
                                queue.append(((nx, ny), dist + 1))

        reachable_cells = len(visited)
        total_cells = n * n
        
        # Check if the probability condition is met
        # reachable_cells / total_cells == 33 / 50 is the same as
        # reachable_cells * 50 == total_cells * 33
        if reachable_cells * 50 == total_cells * 33:
            probability = reachable_cells / total_cells
            print(f"{reachable_cells} / {total_cells} = {probability}")
            return n

        # Move to the next potential value for n
        n += 10
        # Add a safeguard to prevent an infinite loop
        if n > 200:
            print("Failed to find a solution for n up to 200.")
            return None

# Execute the function to find and print the solution
find_grid_size()