import sys

# Setting a higher recursion limit for safety, as this is a recursive solution.
sys.setrecursionlimit(3000)

class GeneralizedHanoiSolver:
    """
    Solves a generalized Tower of Hanoi puzzle where disks start in arbitrary positions.
    The goal is to move all disks to a single target peg in the minimal number of moves
    using a standard recursive algorithm.
    """
    def __init__(self, initial_pegs_dict, target_peg_index):
        # Deep copy the initial state to avoid modifying the original dict
        self.pegs = {p: list(d) for p, d in initial_pegs_dict.items()}
        all_peg_keys = list(self.pegs.keys())
        self.target_peg = target_peg_index
        
        # Determine the total number of disks from the initial configuration
        all_disks = []
        for disk_list in self.pegs.values():
            all_disks.extend(disk_list)
        if not all_disks:
            self.num_disks = 0
        else:
            self.num_disks = max(all_disks)
        
        self.moves = 0

    def find_disk(self, disk_num):
        """Finds the peg where a specific disk is located."""
        for peg_idx, disks in self.pegs.items():
            if disk_num in disks:
                return peg_idx
        raise ValueError(f"Disk {disk_num} could not be found on any peg.")

    def move_disk(self, from_peg, to_peg):
        """
        Executes a single disk move, validates it, updates the state,
        and prints the move.
        """
        # --- Validation ---
        if not self.pegs.get(from_peg):
            raise RuntimeError(f"Error: Cannot move from an empty peg {from_peg}.")
        
        disk_to_move = self.pegs[from_peg][-1]
        
        if self.pegs.get(to_peg) and self.pegs[to_peg] and self.pegs[to_peg][-1] < disk_to_move:
            raise RuntimeError(f"Error: Cannot place disk {disk_to_move} on top of smaller disk "
                               f"{self.pegs[to_peg][-1]} on peg {to_peg}.")

        # --- Execution ---
        disk = self.pegs[from_peg].pop()
        self.pegs[to_peg].append(disk)
        self.moves += 1
        
        # Print the move
        print(f"Move disk {disk} from Peg {from_peg} to Peg {to_peg}")

    def solve(self):
        """Starts the solving process for all disks towards the target peg."""
        print("This script solves the puzzle by moving all disks to the target peg.")
        print("Each move is printed below:\n")
        self.move_tower_to_destination(self.num_disks, self.target_peg)
        print(f"\nTarget position reached on Peg {self.target_peg}: {self.pegs[self.target_peg]}")
        print(f"Minimal amount of moves: {self.moves}")
        return self.moves

    def move_tower_to_destination(self, n, dest_peg):
        """
        The core recursive function. Ensures that disks 1 through n
        end up on the destination peg `dest_peg`.
        """
        if n == 0:
            return

        current_peg = self.find_disk(n)

        # If the largest disk of the current sub-problem (`n`) is already
        # at the destination peg, we only need to worry about moving the
        # smaller disks (1 to n-1) on top of it.
        if current_peg == dest_peg:
            self.move_tower_to_destination(n - 1, dest_peg)
        else:
            # The disk `n` is not at its destination. We must move it.
            # To do so, we must first move all smaller disks (1 to n-1)
            # to a single auxiliary peg.
            
            # Choose an auxiliary peg that is neither the current peg of disk `n`
            # nor its final destination peg.
            available_pegs = list(self.pegs.keys())
            aux_peg = next(p for p in available_pegs if p != current_peg and p != dest_peg)
            
            # Step 1: Move the tower of n-1 disks to the auxiliary peg.
            self.move_tower_to_destination(n - 1, aux_peg)
            
            # Step 2: Now that smaller disks are out of the way, move disk `n`.
            self.move_disk(current_peg, dest_peg)
            
            # Step 3: Move the tower of n-1 disks from the auxiliary peg
            # to the final destination, on top of disk `n`.
            self.move_tower_to_destination(n - 1, dest_peg)

# --- Main execution block ---
# Problem Configuration
start_position = {
    0: [7, 3, 2],
    1: [1],
    2: [8, 6],
    3: [9, 5, 4],
    4: []
}
target_peg_index = 4

solver = GeneralizedHanoiSolver(start_position, target_peg_index)
total_moves = solver.solve()
# The final answer will be extracted from this output.
final_answer = total_moves
# print(f"\n<<< {final_answer} >>>") # For internal verification