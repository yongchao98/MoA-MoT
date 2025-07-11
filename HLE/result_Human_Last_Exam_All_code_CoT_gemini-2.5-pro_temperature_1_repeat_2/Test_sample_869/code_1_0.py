import math

def calculate_probability():
    """
    Calculates the probability as described in the problem.
    
    The steps are:
    1. Calculate S, the total number of possible distributions.
    2. Calculate F, the number of favorable distributions.
    3. Compute the probability P = F / S.
    """
    
    # Step 1: Calculate S, the total number of ways to distribute the items.
    # This is the multinomial coefficient (25 choose 5, 5, 5, 5, 5).
    # S = 25! / (5!^5)
    fact = math.factorial
    S = fact(25) // (fact(5)**5)

    # Step 2: Calculate F, the number of favorable distributions.
    # There are 5! ways to assign a unique dominant type to each person.
    num_assignments = fact(5)

    # For a fixed assignment (Person i is dominant in Type i), we count the ways.
    # Case A: "Pure" distribution (Person i gets 5 of Type i).
    # This corresponds to one specific distribution matrix C.
    # The number of ways to form the hands for this single matrix is 1.
    ways_pure = 1

    # Case B: "3+1+1" distribution (e.g., Person i gets 3 of Type i and 1 of two other types).
    # There are 12 such distribution matrices for a fixed assignment.
    # The number of ways to form the hands for one such matrix is calculated as:
    # P1 gets {3T1, 1T2, 1T5}: C(5,3)*C(5,1)*C(5,1) = 250
    # P2 gets {1T1, 3T2, 1T3}: C(2,1)*C(4,3)*C(5,1) = 40
    # P3 gets {1T2, 3T3, 1T4}: C(1,1)*C(4,3)*C(5,1) = 20 (typo in thought, was C(1,1) for P3 T2) -> C(1,1) is correct for remaining T2s
    # After P1 and P2 hands are dealt, items left: {1T1, 1T2, 4T3, 5T4, 4T5}
    # P3 gets {1T2, 3T3, 1T4}: C(1,1)*C(4,3)*C(5,1) = 1*4*5 = 20
    # After P3 hand is dealt: {1T1, 0T2, 1T3, 4T4, 4T5}
    # P4 gets {1T3, 3T4, 1T5}: C(1,1)*C(4,3)*C(4,1) = 1*4*4 = 16
    # After P4 hand is dealt: {1T1, 0T2, 0T3, 1T4, 3T5}
    # P5 gets {1T1, 1T4, 3T5}: C(1,1)*C(1,1)*C(3,3) = 1
    # Total ways = 250 * 40 * 20 * 16 * 1 = 3,200,000 which is 20**5.
    
    ways_311_per_matrix = 20**5
    num_311_matrices = 12
    ways_311_total = num_311_matrices * ways_311_per_matrix
    
    # Total favorable ways for a fixed assignment
    F_fixed = ways_pure + ways_311_total
    
    # Total F is F_fixed multiplied by the number of dominant type assignments
    F = num_assignments * F_fixed
    
    # Step 3: Calculate the probability P = F / S
    P = F / S
    
    print(f"Total number of ways to distribute the items (S):")
    print(S)
    print("\nNumber of favorable distributions (F):")
    print(F)
    print("\nThe probability is P = F / S:")
    print(f"P = {F} / {S}")
    print(f"P ≈ {P}")

calculate_probability()