import math
from itertools import chain, combinations

def get_subsets(s):
    """Generates all subsets of a given set."""
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))

def is_interval(subset, universe):
    """
    Checks if a subset forms a consecutive interval within a given ordered universe.
    For example, for universe [0, 2, 4], {0, 2} is an interval, but {0, 4} is not.
    """
    if not subset:
        return True
    subset = sorted(subset)
    min_val = subset[0]
    max_val = subset[-1]
    
    start_index = universe.index(min_val)
    end_index = universe.index(max_val)
    
    # The length of the slice from the universe should match the length of the subset
    return len(subset) == (end_index - start_index + 1)

def is_bijunctive(s):
    """
    Checks if a set S defines a bijunctive relation for n=4.
    This is true if S is a union of an interval of even numbers and an interval of odd numbers.
    """
    universe_even = [0, 2, 4]
    universe_odd = [1, 3]
    
    s_even = {x for x in s if x % 2 == 0}
    s_odd = {x for x in s if x % 2 != 0}
    
    return is_interval(s_even, universe_even) and is_interval(s_odd, universe_odd)

def solve():
    """
    Calculates the number of sets S for which the described CSP is NP-hard.
    """
    base_set = {0, 1, 2, 3, 4}
    total_subsets = 2**len(base_set)
    
    p_time_count = 0
    all_possible_s = get_subsets(base_set)
    
    for s_tuple in all_possible_s:
        s = set(s_tuple)
        if is_bijunctive(s):
            p_time_count += 1
            
    np_hard_count = total_subsets - p_time_count
    
    print(f"Let S be a subset of {{0, 1, 2, 3, 4}}.")
    print(f"The total number of possible sets S is 2^{len(base_set)} = {total_subsets}.")
    print(f"A set S leads to a problem in P if it defines a 'bijunctive' relation.")
    print(f"The number of such P-time (bijunctive) sets found is: {p_time_count}.")
    print(f"The number of NP-hard sets is the total number of sets minus the number of P-time sets.")
    print(f"Equation: {total_subsets} - {p_time_count} = {np_hard_count}")

solve()