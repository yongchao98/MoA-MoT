import itertools

class MatchupFinder:
    """
    This class solves the problem of finding the largest set of compatible tennis matchups.
    """
    def __init__(self, num_players, group_size):
        self.num_players = num_players
        self.group_size = group_size
        self.players = list(range(1, self.num_players + 1))
        
        self.all_matchups = []
        self.compat_graph = {}
        self.best_solution_indices = []

    def _generate_all_matchups(self):
        """Generates all possible matchups (combinations of players)."""
        # Using frozenset for efficient intersection operations
        self.all_matchups = [
            frozenset(m) for m in itertools.combinations(self.players, self.group_size)
        ]

    def _build_compatibility_graph(self):
        """Builds a graph where an edge connects two compatible matchups."""
        num_matchups = len(self.all_matchups)
        self.compat_graph = {i: [] for i in range(num_matchups)}
        for i in range(num_matchups):
            for j in range(i + 1, num_matchups):
                # Two matchups are compatible if they share at most 2 players.
                if len(self.all_matchups[i].intersection(self.all_matchups[j])) <= 2:
                    self.compat_graph[i].append(j)
                    self.compat_graph[j].append(i)
    
    def _find_max_clique_recursive(self, potential_indices, current_clique):
        """
        A recursive backtracking function to find the maximum clique.
        
        :param potential_indices: A list of indices of matchups that could extend the current clique.
        :param current_clique: A list of indices forming the current clique.
        """
        # Pruning: if the current branch can't possibly beat the best solution, stop.
        if len(current_clique) + len(potential_indices) <= len(self.best_solution_indices):
            return

        # Check if we have found a new best solution
        if len(current_clique) > len(self.best_solution_indices):
            self.best_solution_indices = list(current_clique)

        # Iterate through each potential candidate to extend the clique
        for i, candidate_idx in enumerate(potential_indices):
            
            # The new set of potential candidates must be compatible with the current candidate_idx
            # and must come from the remaining part of the current potential_indices list.
            new_potential_indices = [
                p_idx for p_idx in potential_indices[i + 1:] 
                if p_idx in self.compat_graph[candidate_idx]
            ]
            
            self._find_max_clique_recursive(new_potential_indices, current_clique + [candidate_idx])

    def solve(self):
        """Executes the plan to find and print the solution."""
        print("Generating all possible matchups...")
        self._generate_all_matchups()
        
        print("Building compatibility graph...")
        self._build_compatibility_graph()
        
        print("Searching for the largest list of compatible matchups...")
        initial_candidates = list(range(len(self.all_matchups)))
        
        # Heuristic: Sort candidates by degree to find large cliques faster, improving pruning.
        initial_candidates.sort(key=lambda i: len(self.compat_graph[i]), reverse=True)

        self._find_max_clique_recursive(initial_candidates, [])
        
        print("\n--- Solution Found ---")
        print(f"The largest list of matchups has a size of: {len(self.best_solution_indices)}")
        print("One such list of matchups is:")
        for index in self.best_solution_indices:
            # Sort players in each matchup for consistent output
            matchup_players = sorted(list(self.all_matchups[index]))
            print(f"- {matchup_players}")

if __name__ == '__main__':
    # Problem parameters
    NUM_PLAYERS = 11
    GROUP_SIZE = 4 # Doubles match

    finder = MatchupFinder(NUM_PLAYERS, GROUP_SIZE)
    finder.solve()