class HanoiSolver:
    """
    Solves a generalized Tower of Hanoi problem by moving the largest disks
    first. It prints each move and calculates the total minimal moves.
    """

    def __init__(self, initial_state):
        self.pegs = {k: list(v) for k, v in initial_state.items()}
        self.move_count = 0
        self.phase_costs = []
        self._current_phase_cost = 0

    def move_disk(self, src_peg, dest_peg):
        """Moves a single disk, validates the move, and prints the action."""
        if not self.pegs[src_peg]:
            raise IndexError(f"Cannot move from empty peg {src_peg}")

        disk = self.pegs[src_peg].pop()

        if self.pegs[dest_peg] and self.pegs[dest_peg][-1] < disk:
            # Put the disk back for state consistency before erroring
            self.pegs[src_peg].append(disk)
            raise ValueError(f"Illegal move: Cannot place disk {disk} on smaller disk {self.pegs[dest_peg][-1]}")

        self.pegs[dest_peg].append(disk)
        self.move_count += 1
        self._current_phase_cost += 1
        # print(f"Move {self.move_count}: Move disk {disk} from peg {src_peg} to peg {dest_peg}")
        # print(f"State: {self.pegs}")

    def move_tower(self, num_disks, src_peg, dest_peg, aux_peg):
        """Recursively moves a tower of `num_disks`."""
        if num_disks <= 0:
            return
        # Note: In a real multi-peg scenario, one might use more than one aux peg.
        # Here we follow the standard single-aux-peg algorithm logic.
        self.move_tower(num_disks - 1, src_peg, aux_peg, dest_peg)
        self.move_disk(src_peg, dest_peg)
        self.move_tower(num_disks - 1, aux_peg, dest_peg, src_peg)

    def start_phase(self):
        """Resets the cost counter for a new phase."""
        self._current_phase_cost = 0

    def end_phase(self):
        """Stores the cost of the completed phase."""
        self.phase_costs.append(self._current_phase_cost)

    def solve(self):
        """Executes the step-by-step solution for the specific puzzle."""
        print("Solving puzzle...")
        print(f"Initial State: {self.pegs}\n")
        
        # Phase 1: Place disk 9
        self.start_phase()
        # Move disk 1 to free up peg 1
        self.move_disk(1, 0)
        # Move tower [5, 4] from peg 3 to peg 1 using peg 4 as auxiliary
        self.move_tower(2, 3, 1, 4)
        # Move disk 9 to its final place
        self.move_disk(3, 4)
        self.end_phase()
        print(f"Phase 1 (Place 9) complete. State: {self.pegs}")

        # Phase 2: Place disk 8
        self.start_phase()
        # Move disk 6 out of the way
        self.move_disk(2, 3)
        # Move disk 8
        self.move_disk(2, 4)
        self.end_phase()
        print(f"Phase 2 (Place 8) complete. State: {self.pegs}")

        # Phase 3: Place disk 7
        self.start_phase()
        # Move tower [3, 2, 1] off disk 7. Pegs 2 and 3 can be used as dest/aux.
        self.move_tower(3, 0, 2, 3)
        # Move disk 7
        self.move_disk(0, 4)
        self.end_phase()
        print(f"Phase 3 (Place 7) complete. State: {self.pegs}")

        # Phase 4: Place disk 6
        self.start_phase()
        self.move_disk(3, 4)
        self.end_phase()
        print(f"Phase 4 (Place 6) complete. State: {self.pegs}")

        # Phase 5: Place disk 5
        self.start_phase()
        # Move disk 4 out of the way to temporary peg 0
        self.move_disk(1, 0)
        # Move disk 5
        self.move_disk(1, 4)
        self.end_phase()
        print(f"Phase 5 (Place 5) complete. State: {self.pegs}")

        # Phase 6: Place disk 4
        self.start_phase()
        self.move_disk(0, 4)
        self.end_phase()
        print(f"Phase 6 (Place 4) complete. State: {self.pegs}")

        # Phase 7: Place disk 3
        self.start_phase()
        # Move tower [2, 1] to temporary peg 0
        self.move_tower(2, 2, 0, 1)
        # Move disk 3
        self.move_disk(2, 4)
        self.end_phase()
        print(f"Phase 7 (Place 3) complete. State: {self.pegs}")

        # Phase 8: Place tower [2, 1]
        self.start_phase()
        self.move_tower(2, 0, 4, 1)
        self.end_phase()
        print(f"Phase 8 (Place 2, 1) complete. State: {self.pegs}")

        print("\nSolution Complete!")
        print(f"Target achieved in {self.move_count} moves.")

        # Print the final equation as requested
        equation_str = " + ".join(map(str, self.phase_costs))
        print(f"Calculation of moves by phase: {equation_str} = {sum(self.phase_costs)}")


# Initial configuration of the pegs
start_position = {
    0: [7, 3, 2],
    1: [1],
    2: [8, 6],
    3: [9, 5, 4],
    4: []
}

# Create a solver and run the solution
solver = HanoiSolver(start_position)
solver.solve()

final_answer = solver.move_count
print(f"<<<{final_answer}>>>")