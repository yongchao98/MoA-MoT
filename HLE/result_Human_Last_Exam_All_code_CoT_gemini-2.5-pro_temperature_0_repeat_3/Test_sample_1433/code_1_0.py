import math

def solve_hypercube_problem():
    """
    This function solves the entire problem as requested.
    It calculates the expected time and variance for d=14, the expected time for d=15,
    and checks the inequality.
    """

    # Part 1 & 2: Calculate EX_14 and Var_14
    d = 14

    # Let E_k be the expected time to reach state 0 (meeting) from a state
    # with Hamming distance k. We need E_d. We know E_0 = 0.
    # Let delta_k = E_k - E_{k-2}.
    # The recurrence relation for delta_k can be derived as:
    # k(k-1) * delta_k - (d-k)(d-k-1) * delta_{k+2} = d^2
    
    deltas = {}  # Using a dictionary to map k -> delta_k

    # Base case for the recurrence at k=d
    # For k=d, the delta_{d+2} term is 0, so d(d-1) * delta_d = d^2
    deltas[d] = d**2 / (d * (d - 1))

    # Solve for delta_k backwards from k=d-2 down to 2
    for k in range(d - 2, 0, -2):
        # Rearranging the recurrence to solve for delta_k:
        # delta_k = (d^2 + (d-k)(d-k-1) * delta_{k+2}) / (k(k-1))
        term1 = d**2 / (k * (k - 1))
        term2_factor = (d - k) * (d - k - 1) / (k * (k - 1))
        term2 = term2_factor * deltas[k + 2]
        deltas[k] = term1 + term2

    # Calculate E_k values. E_k is the sum of deltas up to that point.
    # E_k = sum_{j=1 to k/2} delta_{2j}
    E_values = {}
    E_values[0] = 0
    current_E = 0
    for k in range(2, d + 2, 2):
        current_E += deltas[k]
        E_values[k] = current_E
    
    EX14 = E_values[d]

    # Now, calculate the variance.
    # Let E_k^(2) be the second moment of the time to reach state 0 from k.
    # Let gamma_k = E_k^(2) - E_{k-2}^(2).
    # The recurrence for gamma_k is:
    # k(k-1)/d^2 * gamma_k - (d-k)(d-k-1)/d^2 * gamma_{k+2} = 1 + 2*E_k

    gammas = {} # Using a dictionary to map k -> gamma_k

    # Base case for the recurrence at k=d
    gammas[d] = (d**2 / (d * (d - 1))) * (1 + 2 * E_values[d])

    # Solve for gamma_k backwards from k=d-2 down to 2
    for k in range(d - 2, 0, -2):
        # Rearranging the recurrence to solve for gamma_k:
        # gamma_k = (d^2 * (1 + 2*E_k) + (d-k)(d-k-1) * gamma_{k+2}) / (k(k-1))
        term1_factor = 1 + 2 * E_values[k]
        term1 = (d**2 / (k * (k - 1))) * term1_factor
        term2_factor = (d - k) * (d - k - 1) / (k * (k - 1))
        term2 = term2_factor * gammas[k + 2]
        gammas[k] = term1 + term2
        
    # Calculate E_d^(2), the second moment for the starting distance d.
    # E_d^(2) = sum_{j=1 to d/2} gamma_{2j}
    E2_d = 0
    for k in range(2, d + 2, 2):
        E2_d += gammas[k]
        
    Var14 = E2_d - EX14**2

    # Part 3: Calculate EX_15
    # For d=15 (odd), the Hamming distance between Alice and Bob is always odd.
    # It can never reach 0. Thus, they never meet.
    EX15_str = "inf"

    # Part 4: Check the inequality
    # Is EX_d <= (d/2) * (d^d / d!) for even d?
    # We check for d=14.
    inequality_rhs = (d / 2.0) * (d**d / math.factorial(d))
    inequality_holds = "yes" if EX14 <= inequality_rhs else "no"

    # Print the results in a clear format
    print("EX_14 = " + str(int(EX14)))
    print("D^2X_14 = " + str(int(Var14)))
    print("EX_15 = " + EX15_str)
    print("Is it true that for even d, EX_d <= d/2 * d^d/d!?")
    print("For d = 14:")
    print("EX_14 is approximately " + str(EX14))
    print("The right side of the inequality is d/2 * d^d/d! which is approximately " + str(inequality_rhs))
    print("Since EX_14 is less than or equal to the right side, the answer is " + inequality_holds + ".")

solve_hypercube_problem()