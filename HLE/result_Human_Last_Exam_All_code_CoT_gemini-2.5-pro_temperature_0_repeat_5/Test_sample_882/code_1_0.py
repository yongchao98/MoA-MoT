import math
from scipy.optimize import brentq

def solve_game_theory_problem():
    """
    Solves for the mixed strategy equilibrium probability p and calculates the final result.
    """
    # Step 1: Define game parameters
    N = 8
    M = 3 * N

    # Step 2: Define the two sides of the equilibrium equation E_D(p) = E_S(p).
    # This equation sets the expected payoff for the discrete strategy (LHS) equal to
    # the expected payoff for the uniform split strategy (RHS).

    # E_D(p): Expected payoff for the discrete strategy.
    def lhs(p):
        # This formula calculates E[1/(k+1)] where k is the number of other discrete
        # players in the chosen race. k ~ Binomial(M-1, p/N).
        # The closed form is (1 - (1-q)^(n+1))/((n+1)q) with n=M-1, q=p/N.
        # This simplifies to N*(1-(1-p/N)**M)/(M*p).
        if p == 0:
            # Limit of LHS as p -> 0 is 1/N
            return 1.0 / N
        return N * (1 - (1 - p / N)**M) / (M * p)

    # E_S(p): Expected payoff for the uniform split strategy.
    def rhs(p):
        # This formula calculates N * E[((N-1)/N)^K / (M-K)], where K is the number
        # of opponents playing the discrete strategy. K ~ Binomial(M-1, p).
        if p == 0: # K=0 with probability 1
            return N * (1.0 / M)
        if p == 1: # K=M-1 with probability 1
            return N * (((N - 1.0) / N)**(M - 1))

        total_expectation = 0
        for k in range(M):  # k iterates from 0 to M-1
            # Binomial probability P(K=k) for K ~ Binomial(M-1, p)
            try:
                binom_prob = math.comb(M - 1, k) * (p**k) * ((1 - p)**(M - 1 - k))
            except ValueError:
                # This can happen if p is very close to 0 or 1, causing (1-p) or p to be 0.
                # The loop handles these cases, but this is for robustness.
                continue

            term_value = ((N - 1.0) / N)**k
            term = binom_prob * term_value / (M - k)
            total_expectation += term
        return N * total_expectation

    # Step 3: Define the equation f(p) = 0 to be solved.
    def equation_to_solve(p):
        return lhs(p) - rhs(p)

    # Step 4: Solve the equation for p numerically.
    # We need to find the root in the interval (0, 1).
    # A quick check shows equation_to_solve(epsilon) > 0 and equation_to_solve(1-epsilon) < 0,
    # so a root is guaranteed to exist in the interval.
    try:
        p_solution = brentq(equation_to_solve, 1e-9, 1 - 1e-9, xtol=1e-12, rtol=1e-12)
    except ValueError:
        print("Error: Could not find a solution for p. The function may not cross zero.")
        return

    # Step 5: Calculate and print the final result as per the problem statement.
    # The final equation is result = floor(10000 * (1 - p)).
    # The numbers in this equation are 10000, 1, and the calculated p.
    one_minus_p = 1 - p_solution
    factor = 10000
    final_value = factor * one_minus_p
    result = math.floor(final_value)
    
    print(result)

if __name__ == '__main__':
    solve_game_theory_problem()