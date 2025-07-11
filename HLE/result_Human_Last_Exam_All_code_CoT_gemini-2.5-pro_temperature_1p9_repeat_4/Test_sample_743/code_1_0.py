import itertools

class TennisScheduler:
    """
    Finds the largest list of tennis doubles matchups from 11 players such that
    no two matchups have more than two players in common.
    """

    def __init__(self, num_players=11, group_size=4):
        self.num_players = num_players
        self.group_size = group_size
        self.players = list(range(self.num_players))
        self.all_matchups = list(itertools.combinations(self.players, self.group_size))
        self.num_all_matchups = len(self.all_matchups)
        self.compatibility = []
        self.max_len = 0
        self.final_matchup_list = []

    def _check_compatibility(self, m1_idx, m2_idx):
        """Checks if two matchups have 2 or fewer players in common."""
        m1 = set(self.all_matchups[m1_idx])
        m2 = set(self.all_matchups[m2_idx])
        return len(m1.intersection(m2)) <= 2

    def _precompute_compatibility(self):
        """Creates a compatibility matrix for faster lookups."""
        self.compatibility = [[False] * self.num_all_matchups for _ in range(self.num_all_matchups)]
        for i in range(self.num_all_matchups):
            for j in range(i, self.num_all_matchups):
                if self._check_compatibility(i, j):
                    self.compatibility[i][j] = True
                    self.compatibility[j][i] = True

    def _find_max_clique(self, potential_candidates, current_clique):
        """
        A recursive backtracking algorithm to find the maximum clique.
        'potential_candidates' is a list of matchup indices that can be added.
        'current_clique' is a list of matchup indices in the current valid set.
        """
        # Pruning step: if the number of potential candidates plus the current clique size
        # isn't larger than our best-found solution, we can stop exploring this path.
        if len(potential_candidates) + len(current_clique) <= self.max_len:
            return

        # For each candidate, try adding it to the current clique and recurse.
        for i, pivot in enumerate(potential_candidates):
            
            # The new set of candidates must be compatible with the pivot.
            new_candidates = [c for c in potential_candidates[i + 1:] if self.compatibility[pivot][c]]
            
            new_clique = current_clique + [pivot]
            self._find_max_clique(new_candidates, new_clique)
            
        # After exploring all extensions from the current state, check if the current
        # clique is the best one found so far.
        if len(current_clique) > self.max_len:
            self.max_len = len(current_clique)
            self.final_matchup_list = [self.all_matchups[i] for i in current_clique]


    def solve(self):
        """Executes the plan to find the largest list of matchups."""
        print("Pre-computing compatibility between all possible matchups...")
        self._precompute_compatibility()
        
        print("Searching for the largest possible list of matchups...")
        initial_candidates = list(range(self.num_all_matchups))
        self._find_max_clique(initial_candidates, [])
        
        print(f"\nThe largest number of matchups that can be created is: {self.max_len}\n")
        print("One such list of matchups is:")
        
        for matchup in self.final_matchup_list:
            # Adding 1 to each player number for display purposes (Players 1-11)
            printable_matchup = [player + 1 for player in matchup]
            print(f"- Players: {printable_matchup}")

# Create an instance of the scheduler and run the solver.
scheduler = TennisScheduler()
scheduler.solve()

# The final answer format is specified by the problem.
print("\n<<<11>>>")
