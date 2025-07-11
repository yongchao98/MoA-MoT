import collections

def solve_for_n():
    """
    This script finds the grid size 'n' based on the given movement rules and probability.
    
    The script operates as follows:
    1.  It defines a function to calculate the number of reachable cells for a given grid size 'n'.
        This function, `count_reachable_cells(n)`, implements the object's movement rules.
    2.  Movement Rules:
        - The object starts at cell c2, which is (3, 2) in a 1-based coordinate system.
        - A "diagonal move" allows the object to move to any cell on the same diagonal. This costs 1 move.
        - A "border move" is possible only if the object is at a border cell. It allows movement to an adjacent cell that is also on the border. This also costs 1 move.
    3.  A Breadth-First Search (BFS) is used to find all cells reachable within a maximum of 3 moves.
    4.  The probability of selecting a reachable cell is given as 66% (0.66). If R is the number
        of reachable cells, the governing equation is R / (n*n) = 0.66.
    5.  Based on an analysis of the probability R / (n^2) = 33/50, n must be a multiple of 10.
    6.  The script iterates through potential values of n (10, 20, 30, ...) and, for each n, 
        calculates the number of reachable cells R to find which n satisfies the probability.
    7.  Once the correct 'n' is found, the script prints the verification equation and the final answer.
    """

    def count_reachable_cells(n):
        """
        Calculates the total number of cells reachable within 3 moves on an n x n grid
        using a Breadth-First Search (BFS).
        """
        start_pos = (3, 2)  # Cell c2
        q = collections.deque([(start_pos, 0)]) # Queue stores (position, moves_taken)
        visited = {start_pos} # Set of reachable cells, including the start

        while q:
            pos, moves = q.popleft()

            if moves >= 3:
                continue

            (x, y) = pos

            # Rule 1: Diagonal Movement
            for dx, dy in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
                for k in range(1, n + 1):
                    nx, ny = x + k * dx, y + k * dy
                    if 1 <= nx <= n and 1 <= ny <= n:
                        if (nx, ny) not in visited:
                            visited.add((nx, ny))
                            q.append(((nx, ny), moves + 1))
                    else:
                        break
            
            # Rule 2: Border Movement
            is_on_border = (x == 1 or x == n or y == 1 or y == n)
            if is_on_border:
                for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                    nx, ny = x + dx, y + dy
                    if 1 <= nx <= n and 1 <= ny <= n:
                        is_neighbor_on_border = (nx == 1 or nx == n or ny == 1 or ny == n)
                        if is_neighbor_on_border:
                            if (nx, ny) not in visited:
                                visited.add((nx, ny))
                                q.append(((nx, ny), moves + 1))

        return len(visited)

    target_probability = 0.66
    
    # Based on probability analysis, n must be an even multiple of 10.
    # We will test n = 10, 20, 30, ...
    for m in range(1, 11):
        n = 10 * m
        total_cells = n * n
        reachable_cells = count_reachable_cells(n)
        probability = reachable_cells / total_cells
        
        # Check if the calculated probability matches the target
        if abs(probability - target_probability) < 0.001:
            print(f"Searching for grid size n where the probability of reaching a cell is {target_probability:.2f}...")
            print(f"Found a match for n = {n}.")
            print("-" * 30)
            print(f"For a grid of size {n}x{n}, the total number of cells is {total_cells}.")
            print(f"The number of cells reachable from c2 within 3 moves is {reachable_cells}.")
            print("The resulting probability is calculated as:")
            print(f"Reachable Cells / Total Cells = {reachable_cells} / {total_cells} = {probability:.2f}")
            print("-" * 30)
            print(f"This matches the required probability of 66%. Therefore, the value of n is {n}.")
            print(f"<<<{n}>>>")
            return
            
    print("Could not find a value of n that satisfies the condition within the test range.")

# Run the solver
solve_for_n()