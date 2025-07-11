import math

def solve_expected_value():
    """
    This function calculates the expected value based on the problem description.
    """
    # Step 1: Determine h.
    # From the analysis of the recursive structure f, the smallest value of x
    # for a guaranteed traversal is 0.5.
    h = 0.5

    # Step 2: Determine probabilities for g and solve for P(g).
    # Probability of the direct path in g.
    p1_g = h
    # Probability of the path leading to six g's.
    # p1 + p2 + p3 = 1 and p2 = 2*p3 => h + 3*p3 = 1 => p3 = (1-h)/3
    p3_g = (1.0 - h) / 3.0
    
    # Solve P(g) = p1_g + p3_g * P(g)^6 numerically.
    # Start with an initial guess and iterate for convergence.
    prob_g = 0.5
    for _ in range(100):
        prob_g = p1_g + p3_g * (prob_g ** 6)

    # Step 3: Determine the probability of traversing k.
    # k is a chain of 4 instances of g.
    prob_k = prob_g ** 4

    # Step 4: Calculate the probability of the opponent winning (you losing).
    # This is a binomial probability problem: n=100 trials, p=prob_k success rate.
    # Opponent wins if successes < 6.
    n = 100
    p = prob_k
    
    prob_you_lose = 0.0
    for k_successes in range(6):
        # Calculate P(X=k_successes) for the binomial distribution
        term = math.comb(n, k_successes) * (p ** k_successes) * ((1.0 - p) ** (n - k_successes))
        prob_you_lose += term

    # Step 5: Calculate the final expected value.
    # EV = P(Win)*($1) + P(Lose)*(-$1)
    prob_you_win = 1.0 - prob_you_lose
    win_amount = 1.0
    loss_amount = -1.0
    
    expected_value = prob_you_win * win_amount + prob_you_lose * loss_amount

    # As requested, print the numbers in the final equation.
    print(f"Probability of success for one traversal of k: {p:.6f}")
    print(f"Probability you win (X >= 6): {prob_you_win:.6f}")
    print(f"Probability you lose (X < 6): {prob_you_lose:.6f}")
    print("\nFinal Equation:")
    print(f"Expected Value = {prob_you_win:.6f} * (${win_amount:.2f}) + {prob_you_lose:.6f} * (${loss_amount:.2f})")
    
    print("\nResult:")
    print(f"Your expected value of playing one game is ${expected_value:.2f}")

solve_expected_value()