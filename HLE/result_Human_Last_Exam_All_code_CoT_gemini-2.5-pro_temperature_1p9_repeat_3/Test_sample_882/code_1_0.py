import math
from scipy.optimize import brentq

def find_equilibrium_p():
    """
    Finds the equilibrium probability 'p' by solving the equilibrium condition numerically.
    The equilibrium condition is that the expected payoff for playing the discrete strategy
    is equal to the expected payoff for playing the spread strategy.
    f(p) = E[U_D] - E[U_S] = 0
    """
    N = 8
    # q is the probability of a discrete player NOT choosing a specific race
    q = 1.0 - 1.0/N
    # Total players are 3N. We consider one player, so there are 3N-1 opponents.
    opponents = 3 * N - 1

    # Define the function f(p) = 0 that we need to solve for p.
    def f(p):
        if p < 0 or p > 1:
            return float('inf') # Constraint p to be a probability

        expected_payoff_diff = 0.0
        # Iterate over all possible numbers of discrete opponents 'j'
        for j in range(opponents + 1):
            
            # Binomial probability of having j discrete opponents
            # P(j) = C(opponents, j) * p^j * (1-p)^(opponents-j)
            try:
                # Use logarithms for numerical stability with small probabilities
                log_binom_prob = math.lgamma(opponents + 1) - math.lgamma(j + 1) - math.lgamma(opponents - j + 1)
                log_p_term = j * math.log(p) if p > 0 else (float('-inf') if j > 0 else 0)
                log_1_p_term = (opponents - j) * math.log(1 - p) if p < 1 else (float('-inf') if opponents - j > 0 else 0)
                binom_prob = math.exp(log_binom_prob + log_p_term + log_1_p_term)
            except (ValueError, OverflowError):
                # Fallback to direct computation for cases log can't handle (e.g., p=0 or p=1)
                binom_prob = math.comb(opponents, j) * (p**j) * ((1-p)**(opponents - j))


            # Payoff if playing Discrete Strategy, given j discrete opponents
            # U_D(j) = (N / (j + 1)) * (1 - q^(j + 1))
            ud_j = (N / (j + 1.0)) * (1.0 - q**(j + 1.0))

            # Payoff if playing Spread Strategy, given j discrete opponents
            # U_S(j) = (N / (3N - j)) * q^j
            us_j = (N / (3.0 * N - j)) * (q**j)
            
            payoff_diff = ud_j - us_j
            expected_payoff_diff += binom_prob * payoff_diff

        return expected_payoff_diff

    # Find the root of f(p) in the interval (0, 1).
    # We know f(0) > 0 and f(1) < 0 for N=8, so a root is guaranteed.
    p_equilibrium, r = brentq(f, 0.0, 1.0, full_output=True)
    if not r.converged:
        print("Warning: Root finding did not converge.")

    return p_equilibrium

# --- Main execution ---
p_eq = find_equilibrium_p()
value = 10000 * (1 - p_eq)
final_answer = math.floor(value)

print(f"For N=8, the equilibrium condition is met at p = {p_eq:.8f}")
print("We need to compute floor(10000 * (1-p))")
print(f"The numbers in the final equation are: 10000, 1, {p_eq:.8f}")
print(f"Calculation: floor(10000 * (1 - {p_eq:.8f})) = floor({value:.8f}) = {final_answer}")
print(f"The final answer is: {final_answer}")
<<<94>>>