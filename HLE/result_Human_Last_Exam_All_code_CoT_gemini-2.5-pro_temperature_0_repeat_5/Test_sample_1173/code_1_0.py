import math

def solve():
    """
    This function formalizes the argument that theta = 1/2.
    The core of the argument is that the probability of the stopping time being less than n,
    P(tau < n), can be bounded by a term proportional to n^(-1/2).
    This leads to E[tau] >= n - c*n^(1/2), making theta = 1/2.
    """
    
    # The problem is to find the largest theta (multiple of 1/8) such that
    # E[tau] >= n - c * n^theta
    
    # From the derivation in the text:
    # 1. E[tau] >= n - n * P(tau < n)
    # 2. P(tau < n) <= P(K'_{n-1} >= n^(1/2) - 1)
    #    where K'_{n-1} ~ Binomial(n-1, q = 1/(2*n^(1/2)))
    # 3. Using Cantelli's inequality on K'_{n-1}, we found:
    #    P(K'_{n-1} >= n^(1/2) - 1) <= C * n^(-1/2) for some constant C.
    # 4. This implies E[tau] >= n - n * (C * n^(-1/2)) = n - C * n^(1/2).
    
    # This shows that the inequality holds for theta = 1/2.
    # The question asks for the largest multiple of 1/8.
    # 1/2 = 4/8.
    
    theta_numerator = 4
    theta_denominator = 8
    theta = theta_numerator / theta_denominator
    
    print("Step 1: We establish the bound E[tau] >= n - n * P(tau < n).")
    print("Step 2: We use stochastic domination to bound P(tau < n) by the tail probability of a binomial variable K'_{n-1}.")
    print("P(tau < n) <= P(K'_{n-1} >= n^(1/2) - 1) where K'_{n-1} ~ Binomial(n-1, 1/(2*n^(1/2))).")
    print("Step 3: We apply Cantelli's (one-sided Chebyshev) inequality to bound this probability.")
    print("The mean of K'_{n-1} is mu ~ n^(1/2)/2.")
    print("The variance of K'_{n-1} is sigma^2 ~ n^(1/2)/2.")
    print("The deviation is a ~ n^(1/2)/2.")
    print("P(K'_{n-1} - mu >= a) <= sigma^2 / (sigma^2 + a^2) ~ (n^(1/2)/2) / (n/4) = 2*n^(-1/2).")
    print("Step 4: Substituting this back, we get E[tau] >= n - n * (2*n^(-1/2)) = n - 2*n^(1/2).")
    print(f"Step 5: This proves the inequality for theta = 1/2.")
    print(f"As a multiple of 1/8, theta = {theta_numerator}/{theta_denominator}.")
    
    # The final answer is the value of theta.
    final_theta = 1/2
    print(f"The largest provable theta is {final_theta}.")

solve()
