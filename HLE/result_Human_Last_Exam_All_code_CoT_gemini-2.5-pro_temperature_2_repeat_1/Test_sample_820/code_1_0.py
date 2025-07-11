import collections

def solve_hanoi_puzzle():
    """
    Solves a generalized Tower of Hanoi puzzle using Breadth-First Search (BFS)
    to find the solution with the minimum number of moves.
    """
    
    # Define the starting and target configurations of the disks on the pegs.
    # The state is represented as a tuple of tuples, where each inner tuple is a peg.
    # This format is hashable, allowing it to be stored in the 'visited' set.
    initial_state = (
        (7, 3, 2),    # Peg 0
        (1,),         # Peg 1
        (8, 6),       # Peg 2
        (9, 5, 4),    # Peg 3
        ()            # Peg 4
    )
    target_state = (
        (),
        (),
        (),
        (),
        (9, 8, 7, 6, 5, 4, 3, 2, 1)
    )

    # Initialize a queue for the BFS algorithm.
    # Each item in the queue is a tuple containing the current state and the number of moves taken to reach it.
    queue = collections.deque([(initial_state, 0)])
    
    # A set to keep track of states that have already been visited to avoid cycles and redundant work.
    visited = {initial_state}
    
    # Process the queue until it is empty or the solution is found.
    while queue:
        current_state_tuple, moves = queue.popleft()

        # Check if the current state matches the target state.
        if current_state_tuple == target_state:
            # If a match is found, we have reached the solution in the minimal number of moves.
            # The instruction "output each number in the final equation" is interpreted as
            # printing the final resulting number of moves.
            print("Final equation for minimal moves:")
            print(f"Total Moves = {moves}")
            return

        # For generating next moves, it's easier to work with a mutable list of lists.
        current_state_list = [list(peg) for peg in current_state_tuple]
        num_pegs = len(current_state_list)

        # Iterate over every possible source peg.
        for source_peg_idx in range(num_pegs):
            # We can only move a disk from a non-empty peg.
            if not current_state_list[source_peg_idx]:
                continue
            
            # The disk to be moved is the one on top of the source peg.
            disk_to_move = current_state_list[source_peg_idx][-1]

            # Iterate over every possible destination peg.
            for dest_peg_idx in range(num_pegs):
                # A disk cannot be moved to its current peg.
                if source_peg_idx == dest_peg_idx:
                    continue

                # A move is legal if the destination peg is empty, or if the top disk on the
                # destination peg is larger than the disk we are moving.
                if not current_state_list[dest_peg_idx] or current_state_list[dest_peg_idx][-1] > disk_to_move:
                    # If the move is legal, create the next state.
                    next_state_list = [list(peg) for peg in current_state_list]
                    
                    # Perform the move by popping the disk from the source and appending to the destination.
                    disk = next_state_list[source_peg_idx].pop()
                    next_state_list[dest_peg_idx].append(disk)

                    # Convert the new state back to the hashable tuple format.
                    next_state_tuple = tuple(tuple(peg) for peg in next_state_list)
                    
                    # If this new state has not been visited yet, add it to the queue and visited set.
                    if next_state_tuple not in visited:
                        visited.add(next_state_tuple)
                        queue.append((next_state_tuple, moves + 1))
    
    # This part is reached only if the queue becomes empty and the target was never found.
    print("No solution was found.")

# Execute the solver function.
solve_hanoi_puzzle()
<<<45>>>