import math
from scipy.optimize import fsolve

def solve_game_theory_problem():
    """
    This script solves the game theory problem by finding the mixed strategy
    equilibrium probability 'p' and then calculating the final required value.
    """
    # Step 1: Define problem parameters
    N = 8
    M = 3 * N

    def E_discrete(p):
        """
        Calculates the expected payoff for the discrete strategy.
        A player wins if they overcome or tie with other discrete players in their chosen race.
        The formula is derived from the expectation E[1/(K+1)] where K is the number
        of opponents in the chosen race, K ~ Binomial(M-1, p/N).
        The result is (1 - (1-p/N)^M) / (3*p).
        """
        if p == 0:
            return 1.0  # The limit as p approaches 0
        return (1 - (1 - p / N)**M) / (3 * p)

    def E_split(k, p):
        """
        Calculates the expected payoff for splitting fuel across k races.
        A player wins if at least one of their chosen races is not selected by
        any player using the discrete strategy. The formula uses the
        Principle of Inclusion-Exclusion.
        """
        M_minus_1 = M - 1
        total_payoff = 0
        for j in range(1, k + 1):
            # Probability that a specific set of j races are all free from discrete players
            prob_j_races_free = (1 - p * j / N)**M_minus_1
            
            term = math.comb(k, j) * prob_j_races_free
            if (j - 1) % 2 == 0:  # j is odd
                total_payoff += term
            else:  # j is even
                total_payoff -= term
        return total_payoff

    # Step 2: Find the optimal deviation k_star by checking payoffs at p=1
    payoffs_at_p1 = {k: E_split(k, 1.0) for k in range(2, N + 1)}
    k_star = max(payoffs_at_p1, key=payoffs_at_p1.get)

    # Step 3: Define the equilibrium equation to solve for p
    def equation_to_solve(p_array):
        """
        The equilibrium condition is when the payoffs are equal:
        E_discrete(p) - E_split(k_star, p) = 0.
        """
        p = p_array[0]
        # Add bounds to guide the numerical solver
        if p <= 0 or p > 1:
            return 1e9
        return E_discrete(p) - E_split(k_star, p)

    # Step 4: Solve for the equilibrium probability p using a numerical solver
    initial_guess = [0.95]  # Start with a guess close to 1
    p_equilibrium = fsolve(equation_to_solve, initial_guess)[0]

    # Step 5: Calculate the final result
    final_answer = math.floor(10000 * (1 - p_equilibrium))

    # Output the details of the solution process and the final answer
    print("The game is defined by N={} races and M={} players.".format(N, M))
    print("The optimal deviation strategy is to split fuel across k*={} races.".format(k_star))
    print("\nThe equilibrium is found by solving for 'p' in the equation: Payoff(Discrete) = Payoff(Split).")
    print("The equation is of the form: (1 - (1 - p/{})**{}) / (3*p) = Sum_{{j=1 to {}}}[(-1)**(j-1) * C({},j) * (1 - p*j/{})**{}]".format(N, M, k_star, k_star, N, M-1))
    
    print("\nSolving this equation numerically gives the equilibrium probability p = {:.6f}".format(p_equilibrium))
    
    # Show the numbers in the final solved equation
    lhs_val = E_discrete(p_equilibrium)
    rhs_val = E_split(k_star, p_equilibrium)
    print("\nAt equilibrium (p={:.6f}), the payoffs are equal:".format(p_equilibrium))
    print("  Payoff(Discrete) = {:.6f}".format(lhs_val))
    print("  Payoff(Split k={}) = {:.6f}".format(k_star, rhs_val))
    
    print("\nFinally, we calculate the required value:")
    print("  p = {}".format(p_equilibrium))
    print("  1 - p = {}".format(1 - p_equilibrium))
    print("  10000 * (1 - p) = {}".format(10000 * (1 - p_equilibrium)))
    print("  floor(10000 * (1 - p)) = {}".format(final_answer))

solve_game_theory_problem()
<<<177>>>