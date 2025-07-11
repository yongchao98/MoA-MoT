import sympy

def solve_markov_chain_transience():
    """
    Calculates the infimum of c for a Markov chain to be transient.

    The method uses Pakes's criterion, which relies on the asymptotic drift
    and variance of the jumps for large states k.
    
    The criterion for transience is alpha > sigma^2 / 2, where:
    - alpha = lim_{k->inf} k * E[X_{n+1} - k | X_n = k]
    - sigma^2 = lim_{k->inf} E[(X_{n+1} - k)^2 | X_n = k]
    """
    
    print("This script finds the infimum of c for which the given Markov chain is transient.")
    print("The transition probabilities for large k are:")
    print("P(k, k-2) = 1/4")
    print("P(k, k+2) = 1/4")
    print("P(k, k-1) = 1/4 - c/k")
    print("P(k, k+1) = 1/4 + c/k")
    print("-" * 30)

    # Step 1: Define jumps and their probabilities
    jumps = [-2, -1, 1, 2]
    # For the calculation of asymptotic values, k -> infinity, so c/k -> 0
    # For drift, we must keep the c/k term. For variance, it vanishes.
    probs_drift = ['1/4', '1/4 - c/k', '1/4 + c/k', '1/4']
    probs_variance = [1/4, 1/4, 1/4, 1/4]

    # Step 2: Calculate alpha
    # mu_k = sum(jump * prob)
    # mu_k = (-2)*(1/4) + (-1)*(1/4 - c/k) + (1)*(1/4 + c/k) + (2)*(1/4)
    # mu_k = -1/2 - 1/4 + c/k + 1/4 + c/k + 1/2 = 2*c/k
    # alpha = lim_{k->inf} k * mu_k = k * (2*c/k) = 2*c
    alpha_coeff = 2
    print("Step 1: Calculate the drift parameter alpha.")
    print("The expected drift for large k is mu_k = 2*c/k.")
    print(f"alpha = lim_{{k->inf}} k * mu_k = {alpha_coeff}*c")
    print("-" * 30)

    # Step 3: Calculate sigma^2
    # sigma^2 = lim_{k->inf} sum(jump^2 * prob)
    # The c/k terms vanish as k -> infinity.
    sigma_sq_calc = (jumps[0]**2 * probs_variance[0] + 
                     jumps[1]**2 * probs_variance[1] + 
                     jumps[2]**2 * probs_variance[2] + 
                     jumps[3]**2 * probs_variance[3])
    
    print("Step 2: Calculate the asymptotic jump variance sigma^2.")
    print(f"sigma^2 = lim_{{k->inf}} [(-2)^2*P(k,k-2) + (-1)^2*P(k,k-1) + (1)^2*P(k,k+1) + (2)^2*P(k,k+2)]")
    print(f"sigma^2 = ({jumps[0]**2})*({probs_variance[0]}) + ({jumps[1]**2})*({probs_variance[1]}) + ({jumps[2]**2})*({probs_variance[2]}) + ({jumps[3]**2})*({probs_variance[3]})")
    sigma_sq_val = float(sigma_sq_calc)
    print(f"sigma^2 = {sigma_sq_val}")
    print("-" * 30)

    # Step 4: Apply the transience criterion: alpha > sigma^2 / 2
    print("Step 3: Apply the transience criterion alpha > sigma^2 / 2 and solve for c.")
    print("The inequality for transience is:")
    print(f"{alpha_coeff} * c > {sigma_sq_val} / 2")

    threshold_rhs = sigma_sq_val / 2
    print(f"Which simplifies to:")
    print(f"{alpha_coeff} * c > {threshold_rhs}")

    infimum_c = threshold_rhs / alpha_coeff
    print("Solving for c, we get:")
    print(f"c > {threshold_rhs} / {alpha_coeff}")
    print(f"c > {infimum_c}")
    print("-" * 30)

    # Step 5: Conclude the infimum
    print("The set of c for which the chain is transient is (0.625, infinity).")
    print(f"The infimum of this set is {infimum_c}.")

solve_markov_chain_transience()
<<<0.625>>>