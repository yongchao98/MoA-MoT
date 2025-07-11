def solve():
    """
    This function solves the mathematical problem based on the properties of Stiefel-Whitney classes.

    1. A manifold is "full" if and only if all its Stiefel-Whitney classes w_k are zero for k > 0.
    2. The Stiefel-Whitney classes of a connected sum M = M1 # M2 # ... # Ml are trivial
       if and only if the Stiefel-Whitney classes of each component Mi are trivial.
    3. Therefore, the connected sum is "full" if and only if each component manifold is "full".
    4. The problem requires finding a set of non-full manifolds whose connected sum is full.
    5. This is a logical contradiction. No such set of manifolds exists.
    6. The problem asks for the lexicographically least tuple representing such a set.
       Since the set of solutions is empty, the answer is an empty tuple.
    """
    # The tuple is empty because the problem's conditions are contradictory.
    result_tuple = ()

    # The problem asks to output the flat tuple.
    # For an empty tuple, this is "()".
    print("()")

solve()