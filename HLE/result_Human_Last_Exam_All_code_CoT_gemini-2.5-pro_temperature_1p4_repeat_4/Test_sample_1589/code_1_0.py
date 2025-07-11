import sys

# It is recommended to increase the recursion limit for deep recursion,
# as the search space can be large. For n=10, this is a good precaution.
sys.setrecursionlimit(2000)

class Counter:
    """A simple mutable integer class to pass by reference."""
    def __init__(self):
        self.value = 0

    def increment(self):
        self.value += 1

    def get_value(self):
        return self.value

def find_a(n):
    """
    Calculates a(n) by initiating a backtracking search.
    
    The function a(n) counts the number of permutations y of {n+1, ..., 2n}
    such that the 2n values y_i + i and y_i - i are distinct for i=1..n.

    Args:
        n: The integer for which to calculate a(n).

    Returns:
        The integer value of a(n).
    """
    # The set of numbers to be paired with {1, ..., n}
    y_domain = list(range(n + 1, 2 * n + 1))
    
    # A container for the count, to be passed by reference through the recursion.
    solution_count = Counter()
    
    # Start the backtracking search from the first number, k=0 (representing i=1).
    backtrack(n, 0, set(), set(), y_domain, solution_count)
    
    return solution_count.get_value()

def backtrack(n, k, used_y_values, calculated_values, y_domain, count_container):
    """
    Performs the recursive backtracking search.
    
    Args:
        n: The target number for the problem.
        k: The current 0-indexed number (k represents i=k+1) we are pairing.
        used_y_values: A set of y-values already used in the current path.
        calculated_values: A set of {y_j + (j+1), y_j - (j+1)} generated so far.
        y_domain: The list of possible values for any y_i, i.e., {n+1, ..., 2n}.
        count_container: The object to increment when a full solution is found.
    """
    # If we have successfully found pairs for all k from 0 to n-1, we've found a solution.
    if k == n:
        count_container.increment()
        return

    # The number from the first set is x = k + 1.
    x = k + 1
    
    # Iterate through all possible y candidates from the y_domain.
    for y_candidate in y_domain:
        # If the candidate is already used in this permutation, skip it.
        if y_candidate in used_y_values:
            continue

        # Calculate the two new values based on the current pair (x, y_candidate).
        val1 = y_candidate + x
        val2 = y_candidate - x

        # Check for conflicts with previously calculated values.
        if val1 not in calculated_values and val2 not in calculated_values:
            # If no conflict, add the current choices and recurse to the next level.
            used_y_values.add(y_candidate)
            calculated_values.add(val1)
            calculated_values.add(val2)
            
            backtrack(n, k + 1, used_y_values, calculated_values, y_domain, count_container)
            
            # Backtrack: undo the choices for this path to explore other possibilities.
            calculated_values.remove(val2)
            calculated_values.remove(val1)
            used_y_values.remove(y_candidate)

# Main part of the script
if __name__ == "__main__":
    n = 10
    result = find_a(n)
    
    # The prompt asks to output each number in the final equation.
    # The equation is a(n) = result.
    print(f"a({n}) = {result}")

<<<48>>>