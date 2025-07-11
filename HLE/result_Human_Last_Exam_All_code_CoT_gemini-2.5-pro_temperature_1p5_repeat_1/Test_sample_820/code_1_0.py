import collections

def solve_hanoi_puzzle():
    """
    Solves the specified Tower of Hanoi variant using Breadth-First Search (BFS)
    to find the minimal number of moves.
    """
    # In each peg's list, the disk at the end of the list is the top disk.
    # Disks must be placed in ascending order of size (larger numbers are larger disks).
    # A smaller disk can be placed on a larger disk.
    start_state = (
        (7, 3, 2),  # Peg 0
        (1,),       # Peg 1
        (8, 6),     # Peg 2
        (9, 5, 4),  # Peg 3
        ()          # Peg 4 (target)
    )

    goal_state = (
        (),
        (),
        (),
        (),
        (9, 8, 7, 6, 5, 4, 3, 2, 1)
    )

    # The queue will store tuples of (current_state, path_list)
    queue = collections.deque([(start_state, [])])
    
    # A set to store visited states to prevent cycles and redundant computations.
    # States must be hashable, which is why we use tuples of tuples.
    visited = {start_state}
    
    num_pegs = len(start_state)

    while queue:
        current_state, path = queue.popleft()

        # Check if the current state is the goal state
        if current_state == goal_state:
            num_moves = len(path)
            # Create the equation string "1 + 1 + ... = total"
            equation_str = " + ".join(['1'] * num_moves)
            print(f"{equation_str} = {num_moves}")
            return num_moves

        # Generate all possible next states from the current state
        for source_peg_idx in range(num_pegs):
            # If the source peg is empty, we can't move anything from it
            if not current_state[source_peg_idx]:
                continue
            
            # The top disk is the last element in the tuple
            disk_to_move = current_state[source_peg_idx][-1]

            for dest_peg_idx in range(num_pegs):
                # A disk cannot be moved to the same peg
                if source_peg_idx == dest_peg_idx:
                    continue

                # Check if the move is legal:
                # 1. The destination peg is empty.
                # 2. Or, the disk to move is smaller than the top disk on the destination peg.
                if not current_state[dest_peg_idx] or disk_to_move < current_state[dest_peg_idx][-1]:
                    
                    # Create the new state after the move
                    new_state_list = [list(peg) for peg in current_state]
                    
                    # Perform the move
                    moved_disk = new_state_list[source_peg_idx].pop()
                    new_state_list[dest_peg_idx].append(moved_disk)
                    
                    # Convert the list of lists back to a tuple of tuples to store in the 'visited' set
                    new_state = tuple(tuple(peg) for peg in new_state_list)

                    if new_state not in visited:
                        visited.add(new_state)
                        # A move is defined as (disk, from, to) for tracking purposes
                        new_path = path + [(disk_to_move, source_peg_idx, dest_peg_idx)]
                        queue.append((new_state, new_path))
                        
    # This part should not be reached if a solution exists
    print("No solution found.")

if __name__ == '__main__':
    solve_hanoi_puzzle()
<<<47>>>