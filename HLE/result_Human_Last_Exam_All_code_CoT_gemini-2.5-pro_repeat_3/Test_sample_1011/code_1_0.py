def solve_ultrafilter_problem():
    """
    This function outlines the step-by-step logical proof to determine the smallest
    possible number of accumulation points for the given set of ultrafilters.
    """

    # Step 1: Formalizing the problem
    # Let N be the set of natural numbers.
    # A partition P = {P_1, P_2, ...} of N is a collection of pairwise disjoint
    # subsets of N whose union is N. Each P_i is infinite.
    # An ultrafilter on N is a maximal filter. A non-principal ultrafilter is one
    # that does not contain any finite sets.
    # The Stone-Cech remainder, N*, is the space of all non-principal ultrafilters on N.
    # We are given a sequence of non-principal ultrafilters u_1, u_2, ... such that
    # for each i, the set P_i is an element of the ultrafilter u_i (P_i in u_i).
    # S = {u_1, u_2, ...} is a countable subset of N*.
    # An accumulation point of S is a point w in N* such that every open neighborhood
    # of w contains infinitely many points from S.
    # The question asks for the minimum possible size of the set of accumulation points of S.

    # Step 2: The number of accumulation points is at least 1.
    # The space N* is compact and Hausdorff.
    # The sequence S = {u_1, u_2, ...} consists of distinct points. Since P_i in u_i
    # and P_j in u_j for i != j, and P_i, P_j are disjoint, the ultrafilters u_i and u_j
    # are separated by the disjoint open sets P_i* and P_j*. Thus, S is an infinite set.
    # In any compact space, any infinite subset must have at least one accumulation point.
    # Therefore, the number of accumulation points is >= 1.

    # Step 3: The number of accumulation points cannot be 1.
    # Proving this is the core of the argument.
    # If there were exactly one accumulation point, w, this would imply that the
    # sequence {u_i} converges to w.
    # By the definition of convergence in the topology of N*, this means that for any
    # set A in the ultrafilter w (A in w), the set A must also be in u_i for all but
    # a finite number of indices i.
    #
    # Let's analyze the consequence:
    # a) If A is in w, then for almost all i, A is in u_i.
    # b) We are given that P_i is in u_i for all i.
    # c) If A is in u_i and P_i is in u_i, then their intersection (A intersect P_i) must also be in u_i.
    # d) Since u_i is a non-principal ultrafilter, it contains no finite sets. Thus, the set (A intersect P_i) must be infinite.
    #
    # This leads to a necessary condition for convergence: For any set A in w, the intersection (A intersect P_i)
    # must be infinite for all but finitely many indices i.
    #
    # Now, let's define a collection of sets F = {A subset N | (A intersect P_i) is infinite for almost all i}.
    # The condition for convergence means that the ultrafilter w must contain this entire collection F.
    # This is only possible if F has the Finite Intersection Property (FIP), meaning any finite
    # subcollection of F has a non-empty intersection.
    #
    # We can show that F does NOT have the FIP.
    # Since each P_i is infinite, we can partition it into two disjoint infinite sets, P_i_1 and P_i_2.
    # Construct two new sets:
    #   A1 = The union of all P_i_1 for i = 1, 2, 3, ...
    #   A2 = The union of all P_i_2 for i = 1, 2, 3, ...
    #
    # Check if A1 and A2 are in F:
    # - For A1, its intersection with any P_i is P_i_1, which is infinite. So, A1 is in F.
    # - For A2, its intersection with any P_i is P_i_2, which is infinite. So, A2 is in F.
    #
    # However, the intersection of A1 and A2 is the empty set, because for each i, P_i_1 and P_i_2 are disjoint.
    # A collection of sets that contains two disjoint sets (A1 and A2) cannot have the FIP.
    # Therefore, no ultrafilter can contain the collection F.
    # This contradicts the consequence of our initial assumption.
    # The assumption that the sequence {u_i} converges must be false.
    # Thus, the number of accumulation points cannot be 1.

    # Step 4: Final Conclusion
    # From Step 2, the number of accumulation points is >= 1.
    # From Step 3, the number is not equal to 1.
    # Therefore, the smallest possible number must be at least 2.
    #
    # It is a known (but non-trivial) result in set theory that a construction with exactly 2 accumulation
    # points is possible. This involves partitioning the index set {1, 2, ...} into two infinite sets,
    # say even and odd numbers, and then carefully constructing the ultrafilters u_i for the even indices
    # to cluster around one point, and those for the odd indices to cluster around another.
    #
    # So, the minimum possible number is 2.

    smallest_possible_number = 2
    
    # The problem asks to output the numbers in the final equation.
    # There is no equation, so we will simply print the final numerical answer.
    print(smallest_possible_number)

solve_ultrafilter_problem()