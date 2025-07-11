import collections

def solve():
    """
    Finds the value of n by searching through even integers and calculating
    the number of reachable cells using a Breadth-First Search (BFS).
    """

    start_n = 4  # The start cell c2 (3,2) must be in the grid.

    for n in range(start_n, 101, 2):  # Iterate through even numbers
        start_pos = (3, 2)
        # A queue for BFS: stores (position, moves_taken)
        q = collections.deque([(start_pos, 0)])
        # A set to keep track of visited cells to avoid cycles and redundant work
        visited = {start_pos}

        # The BFS explores the grid up to 3 moves
        queue_snapshot = list(q)
        while queue_snapshot:
            pos, moves = queue_snapshot.pop(0)

            if moves >= 3:
                continue

            px, py = pos

            # 1. Diagonal Moves
            # Explore all 4 diagonal directions
            for dx, dy in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
                nx, ny = px, py
                while True:
                    nx, ny = nx + dx, ny + dy
                    # Stop if the move goes out of bounds
                    if not (1 <= nx <= n and 1 <= ny <= n):
                        break
                    
                    if (nx, ny) not in visited:
                        visited.add((nx, ny))
                        q.append(((nx, ny), moves + 1))
            
            # 2. Border Moves
            # A border move is only possible if the current cell is on the border
            if px == 1 or px == n or py == 1 or py == n:
                # Explore all 4 adjacent directions (up, down, left, right)
                for bdx, bdy in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
                    nx, ny = px + bdx, py + bdy
                    # Check if the adjacent cell is valid and also on the border
                    if (1 <= nx <= n and 1 <= ny <= n) and (nx == 1 or nx == n or ny == 1 or ny == n):
                        if (nx, ny) not in visited:
                            visited.add((nx, ny))
                            q.append(((nx, ny), moves + 1))
            
            # After processing one level of moves, take a new snapshot of the queue
            if not queue_snapshot:
                 queue_snapshot = [item for item in q if item not in visited]


        num_reachable = len(visited)
        total_cells = n * n
        probability = num_reachable / total_cells

        # Check if the calculated probability matches the given 66%
        if abs(probability - 0.66) < 1e-9: # Using a small tolerance for float comparison
            print(f"For n = {n}:")
            print(f"The total number of cells in the grid is {n} x {n} = {total_cells}.")
            print(f"The number of cells reachable by the object within 3 moves is {num_reachable}.")
            print("The probability of selecting a reachable cell is:")
            print(f"{num_reachable} / {total_cells} = {probability}")
            print(f"\nThus, the value of n is {n}.")
            
            # To meet the output format request
            final_n = n
            final_reachable = num_reachable
            final_total = total_cells
            final_prob_num = int(probability * 100)
            final_prob_den = 100

            print("\nFinal Equation:")
            print(f"{final_reachable} / {final_total} = {final_prob_num} / {final_prob_den}")
            return

solve()
<<<50>>>