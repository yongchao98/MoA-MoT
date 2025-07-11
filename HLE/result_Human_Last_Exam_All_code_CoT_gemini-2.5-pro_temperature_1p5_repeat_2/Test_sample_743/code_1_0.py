import itertools

class TennisScheduler:
    """
    This class solves the tennis matchup problem by finding the largest set of
    4-player groups from 11 players, where no two groups share more than 2 players.
    """

    def __init__(self, num_players=11, group_size=4):
        self.num_players = num_players
        self.group_size = group_size
        self.players = list(range(self.num_players))
        self.all_matchups = [frozenset(c) for c in itertools.combinations(self.players, self.group_size)]
        self.best_solution = []

    def _find_max_sets_recursive(self, current_solution, start_index):
        """
        Recursively search for the largest valid set of matchups using backtracking.
        """
        # Pruning: If the current path can't possibly beat the best solution found so far,
        # stop exploring it.
        if len(current_solution) + (len(self.all_matchups) - start_index) <= len(self.best_solution):
            return

        # Check if the current solution is better than the best one found so far.
        if len(current_solution) > len(self.best_solution):
            self.best_solution = list(current_solution) # Store a copy
            print(f"Found a new largest list with {len(self.best_solution)} matchups.")

        # Iterate through the remaining candidates to extend the current solution.
        for i in range(start_index, len(self.all_matchups)):
            new_candidate = self.all_matchups[i]
            
            # Check for compatibility with all matchups in the current solution.
            is_compatible = True
            for existing_matchup in current_solution:
                if len(new_candidate.intersection(existing_matchup)) > 2:
                    is_compatible = False
                    break
            
            # If compatible, add it and recurse.
            if is_compatible:
                current_solution.append(new_candidate)
                self._find_max_sets_recursive(current_solution, i + 1)
                current_solution.pop() # Backtrack to explore other possibilities.

    def solve(self):
        """
        Executes the search and prints the final result.
        """
        print(f"There are {len(self.all_matchups)} possible matchups in total.")
        print("Searching for the largest valid set of matchups...")
        
        self._find_max_sets_recursive([], 0)

        print("\n--- Search Complete ---")
        if not self.best_solution:
            print("No valid list of matchups could be created.")
        else:
            print(f"\nThe largest list of matchups that can be created contains {len(self.best_solution)} groups.")
            print("\nThe matchups in this list are (players are numbered 0 to 10):")
            # Sort the final list for consistent output
            sorted_solution = sorted([sorted(list(m)) for m in self.best_solution])
            for i, matchup in enumerate(sorted_solution):
                # This line prints each number in the final equation as requested.
                print(f"Match {i+1:2d}: {matchup}")
        
        return len(self.best_solution)

if __name__ == '__main__':
    scheduler = TennisScheduler()
    final_answer = scheduler.solve()
    print(f"\n<<< {final_answer} >>>")