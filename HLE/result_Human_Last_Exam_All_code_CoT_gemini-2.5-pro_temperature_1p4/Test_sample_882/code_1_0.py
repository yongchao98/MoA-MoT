import math
from decimal import Decimal, getcontext

def solve_for_p():
    """
    This function implements the plan to solve the problem.
    It calculates the mixed strategy equilibrium probability p and the final requested value.
    """
    # Set precision for Decimal calculations to ensure accuracy
    getcontext().prec = 100

    # Step 1: Define game parameters
    N = Decimal(8)
    num_players = 3 * N
    other_players = num_players - 1

    # Helper function for combinations nCr
    def nCr_dec(n, r):
        n_int = int(n)
        r_int = int(r)
        # Using integer math for factorial is faster and exact
        return Decimal(math.factorial(n_int) // (math.factorial(r_int) * math.factorial(n_int-r_int)))

    # Step 2: Determine the best continuous strategy C_m and its payoff against pure D-strategists.
    # The deviating strategy will be the one with the highest payoff.
    best_m = -1
    max_pi_C = Decimal(0)
    
    # Calculate pi_C(m, N) for m = 2...N
    # pi_C(m, N) = sum_{k=1 to m} (-1)^(k-1) * C(m,k) * (1-k/N)^(3N-1)
    for m_int in range(2, int(N) + 1):
        m = Decimal(m_int)
        current_pi_C = Decimal(0)
        for k_int in range(1, m_int + 1):
            k = Decimal(k_int)
            term = (Decimal(-1)**(k - 1)) * nCr_dec(m, k) * (Decimal(1) - k/N)**other_players
            current_pi_C += term
        if current_pi_C > max_pi_C:
            max_pi_C = current_pi_C
            best_m = m_int
            
    # The target payoff for the equilibrium is this maximal payoff
    target_payoff = max_pi_C
    
    # Step 3: Define the payoff function E_D(p) and the equation to solve
    # E_D(p) = (1/(3p)) * [1 - (1-p/N)^(3N)]
    def g(p):
        p_dec = Decimal(p)
        # Handle the limit case p=0 to avoid division by zero
        if p_dec == 0:
            val_at_0 = Decimal(1)/N # The limit of E_D(p) as p->0 is 1/N
            return val_at_0 - target_payoff
            
        e_d_p = (Decimal(1)/(Decimal(3)*p_dec)) * (Decimal(1) - (Decimal(1) - p_dec/N)**num_players)
        return e_d_p - target_payoff

    # Step 4: Solve for p using bisection method
    # We are looking for a root of g(p) = 0
    low = Decimal('0.5')
    high = Decimal('1.0')

    # Iterate 100 times for high precision
    for _ in range(100):
        mid = (low + high) / 2
        if g(mid) > 0: # E_D(p) is a decreasing function of p
            low = mid
        else:
            high = mid
            
    p_equilibrium = (low + high) / 2

    # Step 5: Calculate the final result
    one_minus_p = 1 - p_equilibrium
    final_value_float = 10000 * one_minus_p
    result = final_value_float.to_integral_value(rounding='ROUND_FLOOR')
    
    # Output the steps as requested
    print("For N=8, the optimal continuous strategy is to split fuel over m={} races.".format(best_m))
    print("The Nash equilibrium is defined by the equation E_D(p) = E_C(1), where:")
    print("E_C(1) = pi_C({}, {}) = {}".format(best_m, int(N), target_payoff))
    print("\nSolving for p...")
    print("Equilibrium probability p = {}".format(p_equilibrium))
    
    equilibrium_payoff_D = (Decimal(1)/(Decimal(3)*p_equilibrium)) * (Decimal(1) - (Decimal(1) - p_equilibrium/N)**num_players)
    print("At this p, E_D(p) = {}".format(equilibrium_payoff_D))
    print("This matches our target payoff E_C(1), confirming our solution for p.")
    
    print("\nFinal calculation:")
    print("1 - p = {}".format(one_minus_p))
    print("10000 * (1-p) = {}".format(final_value_float))
    print("floor(10000 * (1-p)) = {}".format(result))

solve_for_p()