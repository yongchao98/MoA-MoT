import sys

def solve_set_problem():
    """
    Solves the problem by analyzing the mathematical properties of the set Σ.

    The problem defines Σ based on the condition A + A ⊆ A × A.
    - A + A is a set of integers (sums).
    - A × A is a set of ordered pairs (from the Cartesian product).

    In standard mathematics, an integer cannot be an element of a set of ordered pairs.
    Therefore, the intersection of (A + A) and (A × A) is always the empty set.
    For the subset condition A + A ⊆ A × A to hold, A + A must be empty.

    The sumset A + A is empty if and only if A is the empty set (A = ∅).

    So, the only set satisfying the condition is A = ∅.

    The set Σ is defined by taking all sets that satisfy the condition and removing ∅ and {2}.
    Let S be the set of all valid A. So, S = {∅}.
    Σ = S \ {∅, {2}} = {∅} \ {∅, {2}} = ∅.

    The set Σ is empty. The problem asks to return 0 in this case.
    """
    
    # The set Σ is empty based on the analysis above.
    is_empty = True
    
    if is_empty:
        result = 0
    else:
        # This part of the code would compute min(max(A)) for A in Σ,
        # but it is unreachable because Σ is empty.
        pass

    print("Step 1: The condition is A + A ⊆ A × A.")
    print("Step 2: A + A contains integers, while A × A contains ordered pairs.")
    print("Step 3: An integer cannot be an ordered pair, so for the condition to hold, A + A must be empty.")
    print("Step 4: A + A is empty if and only if A is empty (A = ∅).")
    print("Step 5: The set of solutions is {∅}. Σ is defined by removing ∅ and {2}.")
    print("Step 6: Σ = {∅} \\ {∅, {2}} = ∅.")
    print("Step 7: Since Σ is empty, the required value is 0.")
    print("\nFinal Result:")
    print(result)

solve_set_problem()