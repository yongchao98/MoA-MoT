import itertools

def solve():
    """
    This script calculates the number of NP-hard variations of a given
    constraint satisfaction problem.

    The problem is defined by a set S, a subset of {0, 1, 2, 3, 4}. A clause
    with 4 boolean inputs is satisfied if the number of TRUE inputs is in S.
    We classify the complexity of this problem for all possible 2^5 = 32 sets S.

    According to Schaefer's Dichotomy Theorem for CSPs where constants are
    available, the problem is in P if the relation is entirely affine, Horn,
    dual-Horn, or bijunctive. Otherwise, it is NP-hard.
    """
    
    # Define the sets S that result in tractable (P-time) problems.
    # Sets are represented as frozensets for use in other sets.
    
    # 1. Affine relations
    affine_sets = {
        frozenset(),                # The empty relation
        frozenset({0, 1, 2, 3, 4}), # The trivial relation (always true)
        frozenset({0}),             # all-false
        frozenset({4}),             # all-true
        frozenset({0, 4}),          # all-equal
        frozenset({1, 3}),          # odd number of true inputs (parity)
        frozenset({0, 2, 4})        # even number of true inputs (negated parity)
    }

    # 2. Horn relations (of the form w(x) <= k) that are not affine
    # The base Horn sets are of form {0,1,...,k}. Those for k=0 and k=4 are already affine.
    horn_only_sets = {
        frozenset({0, 1}),          # w(x) <= 1
        frozenset({0, 1, 2}),       # w(x) <= 2
        frozenset({0, 1, 2, 3})     # w(x) <= 3 (or not all true)
    }

    # 3. Dual-Horn relations (of the form w(x) >= k) that are not affine
    # The base dual-Horn sets are of form {k,...,4}. Those for k=4 and k=0 are affine.
    dual_horn_only_sets = {
        frozenset({1, 2, 3, 4}),    # w(x) >= 1 (or not all false)
        frozenset({2, 3, 4}),       # w(x) >= 2
        frozenset({3, 4})           # w(x) >= 3
    }
    
    # 4. Bijunctive relations: For arity 4, these are w(x)<=1 and w(x)>=3.
    # These correspond to S={0,1} and S={3,4}, which are already counted in the
    # Horn and dual-Horn categories.
    
    # The set of all tractable sets is the union of the above collections.
    tractable_sets = affine_sets.union(horn_only_sets).union(dual_horn_only_sets)

    num_affine = len(affine_sets)
    num_horn_only = len(horn_only_sets)
    num_dual_horn_only = len(dual_horn_only_sets)
    num_tractable = len(tractable_sets)
    
    total_sets = 2**5
    
    num_np_hard = total_sets - num_tractable
    
    print("Step-by-step calculation:")
    print(f"1. Total number of possible sets S is 2^5 = {total_sets}.")
    print(f"2. Number of tractable sets S identified:")
    print(f"   - Affine sets: {num_affine}")
    print(f"   - Horn sets (but not affine): {num_horn_only}")
    print(f"   - Dual-Horn sets (but not affine): {num_dual_horn_only}")
    print(f"   - Bijunctive sets are already included in the above categories.")
    print(f"3. The total number of unique tractable sets is the sum of these counts, as they are disjoint:")
    print(f"   Total tractable = {num_affine} + {num_horn_only} + {num_dual_horn_only} = {num_tractable}")
    print(f"4. The number of NP-hard sets is the total number of sets minus the number of tractable sets:")
    print(f"   Total NP-hard = {total_sets} - {num_tractable} = {num_np_hard}")
    
if __name__ == '__main__':
    solve()
<<<19>>>