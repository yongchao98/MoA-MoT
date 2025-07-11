import math
from itertools import chain, combinations

def get_number_of_digits(n):
    """Calculates the number of decimal digits in a positive integer."""
    if n == 0:
        return 1
    # This works for n > 0, which is the case for this problem
    return math.floor(math.log10(n)) + 1

def concat_numbers(c, d):
    """
    Performs the concatenation of two integers.
    This corresponds to the interpretation of the 'A x A' operation.
    """
    return c * (10 ** get_number_of_digits(d)) + d

def check_set_property(A):
    """
    Checks if a set A satisfies the condition A+A is a subset of C(A),
    where C(A) is the set of numbers formed by concatenating elements of A.
    """
    if not A:
        # The empty set trivially satisfies the condition but is excluded from Sigma
        return True

    # 1. Compute the sum set A + A
    sum_set = {a1 + a2 for a1 in A for a2 in A}

    # 2. Compute the concatenation set C(A)
    concatenation_set = {concat_numbers(c, d) for c in A for d in A}

    # 3. Check if sum_set is a subset of concatenation_set
    return sum_set.issubset(concatenation_set)

def find_solution():
    """
    Searches for a set A in Sigma and returns the minimum max element.
    The logical proof shows that Sigma is empty, so this search will not
    find a solution. The function will return 0 as per the problem.
    """
    # We perform a bounded search to demonstrate the process.
    # The search iterates through M = max(A) from 1 upwards.
    search_limit_M = 20  # A reasonable limit for the search demonstration

    for M in range(1, search_limit_M + 1):
        # For a given M, the set A must contain M. Other elements of A
        # must be chosen from {1, 2, ..., M-1}.
        potential_other_elements = list(range(1, M))
        
        # Iterate through all subsets of the potential other elements
        # This will form the part of A that is smaller than M.
        for r in range(len(potential_other_elements) + 1):
            for subset in combinations(potential_other_elements, r):
                A = set(subset)
                A.add(M)

                # The problem statement defines Sigma by excluding O and {2}
                if A == set() or A == {2}:
                    continue
                
                # Check if this set A has the required property
                if check_set_property(A):
                    # If a set is found, M is its maximum element.
                    # Since we are iterating M upwards, this will be the
                    # minimum possible maximum element. Our proof shows this
                    # part of the code is unreachable.
                    return M
    
    # If the search completes without finding any set, it means Sigma is empty
    # within our search bounds. The mathematical proof confirms it's truly empty.
    return 0

# The logical deduction shows that no such set exists and Sigma is empty.
# The problem specifies to return 0 in this case.
result = find_solution()
print(result)