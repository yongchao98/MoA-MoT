import math

def solve_proportionality_problem():
    """
    Solves for s1 (PJR) and s2 (EJR) based on the problem description.
    """
    
    k = 100  # Committee size
    
    # --- Part 1: Solving for s1 (Proportional Justified Representation) ---
    print("--- Solving for s1 (Proportional Justified Representation) ---")
    print(f"Let n be the total number of voters and k = {k} be the committee size.")
    
    A1 = {'a', 'b', 'c', 'x'}
    A2 = {'a', 'b', 'c', 'y'}
    A3 = {'a', 'b', 'c', 'y'}
    A4 = {'a', 'b', 'c', 'z'}
    A5 = {'a', 'b', 'c', 'z'}
    A6 = {'a', 'b', 'c', 'z'}
    initial_ballots = [A1, A2, A3, A4, A5, A6]
    
    print("\nTo leave voter 1 unsatisfied, the committee W must not contain any candidates from their ballot A(1).")
    print(f"So, W must not contain any of {A1}.")

    # Consider the group of the first 6 voters
    N0_size = len(initial_ballots)
    C0 = set.intersection(*initial_ballots)
    
    print(f"\nConsider the group of the first {N0_size} voters.")
    print(f"The intersection of their approved candidates is C0 = {C0}, with size |C0| = {len(C0)}.")
    
    print("\nThe definition of PJR states that for any group of voters N' of size at least n/k who have a common set of approved candidates C, the committee W must contain at least one candidate from C.")
    print(f"Formally: if |N'| >= n/k and C = intersection(A_i for i in N') is not empty, then C intersect W is not empty.")

    print("\nFor our group of {N0_size} voters, their common candidates are {C0}.")
    print("Since voter 1 must be unsatisfied, W cannot contain 'a', 'b', or 'c'. So, C0 intersect W is empty.")
    print("To avoid violating PJR, the PJR condition must not apply to this group. This means the group size must be less than the n/k threshold.")
    
    print(f"\nTherefore, we must have the inequality: |N'| < n/k")
    print(f"Substituting the values: {N0_size} < n / {k}")
    
    # s1 is the smallest integer n > N0_size * k
    s1 = N0_size * k + 1
    
    print(f"This simplifies to n > {N0_size * k}.")
    print(f"The smallest integer n that satisfies this condition is {s1}.")
    print("This size is achievable by adding n-6 voters with empty ballots and choosing a committee W of 100 new candidates.")
    print(f"So, s1 = {s1}")

    # --- Part 2: Solving for s2 (Extended Justified Representation) ---
    print("\n\n--- Solving for s2 (Extended Justified Representation) ---")
    print("EJR states that for any j in {1,...,k}, it's not the case that there is a group N' of size at least j*n/k who are cohesive on j candidates (|intersection(A_i)| >= j) and are all unsatisfied.")
    
    print("\nWe can construct a scenario where only voter 1 is unsatisfied. We place 'y' and 'z' in the committee W.")
    print("This satisfies voters 2 through 6. Other new voters can be satisfied by other candidates in W.")
    print("So, the set of unsatisfied voters is U = {voter 1}.")
    
    print("\nNow we check the EJR condition for the group N' = {voter 1}. |N'| = 1.")
    A1_size = len(A1)
    print(f"The ballot for voter 1 is A(1) = {A1}, so the size of the intersection is |A(1)| = {A1_size}.")
    
    print("\nFor EJR to hold, for any j, we must have |N'| < j*n/k OR |intersection(A_i)| < j.")
    print(f"Substituting the values for N'={{1}}: 1 < j * n / {k}  OR  {A1_size} < j.")
    
    print("\nThis must hold for all j from 1 to 100. The condition (4 < j) is false for j <= 4, so for these j's, the first part must be true.")
    
    # We need 1 < j*n/k for j=1,2,3,4. The strictest case is j=1.
    j = 1
    print(f"For j = {j}: 1 < {j} * n / {k}  =>  n > {k / j}")
    j = 2
    print(f"For j = {j}: 1 < {j} * n / {k}  =>  n > {k / j}")
    j = 3
    print(f"For j = {j}: 1 < {j} * n / {k}  =>  n > {math.ceil(k / j * 100)/100}") # show some decimals
    j = 4
    print(f"For j = {j}: 1 < {j} * n / {k}  =>  n > {k / j}")
    
    strictest_n_lower_bound = k / 1
    s2 = int(strictest_n_lower_bound) + 1
    
    print(f"\nThe strictest condition comes from j=1, which requires n > {strictest_n_lower_bound}.")
    print(f"The smallest integer n that satisfies this is {s2}.")
    print(f"So, s2 = {s2}")
    
    # --- Final Result ---
    print("\n\n--- Final Answer ---")
    final_answer = (s1, s2)
    print(f"The pair (s1, s2) is: {final_answer}")
    
    return final_answer

if __name__ == '__main__':
    solve_proportionality_problem()