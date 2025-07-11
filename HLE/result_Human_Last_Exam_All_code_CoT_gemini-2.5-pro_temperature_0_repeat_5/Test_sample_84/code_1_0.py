import math

def solve():
    """
    This function determines the value of alpha based on the reasoning provided.
    The reasoning involves transforming the problem domain and applying principles from approximation theory.
    
    Let the two sets of points be S1 = {1, ..., n^2} and S2 = {n^2+1, ..., n^10}.
    Let p_n(x) be a polynomial of degree d_n.
    p_n(i) is in [0, 1] for i in S1.
    p_n(i) is in [2, 3] for i in S2.

    1. Transformation: Let x = y^2. Let q(y) = p_n(y^2).
       For p_n(x) to be a polynomial, q(y) must be an even polynomial.
       The degree of p_n is d_n = deg(q)/2.
       The conditions become:
       - q(y) in [0, 1] for y in {sqrt(1), ..., sqrt(n^2)=n}. Call this interval I1 = [1, n].
       - q(y) in [2, 3] for y in {sqrt(n^2+1), ..., sqrt(n^10)=n^5}. Call this interval I2 = [n, n^5].

    2. Normalization: We need to find the degree of a polynomial that can transition from a small value on I1 to a large value on I2.
       Let's map the combined interval [1, n^5] to [0, 1] with a variable u.
       u = (y - 1) / (n^5 - 1).
       The interval I1, [1, n], maps to J1 = [0, (n-1)/(n^5-1)] which is approximately [0, n^-4].
       The interval I2, [n, n^5], maps to J2 = [(n-1)/(n^5-1), 1] which is approximately [n^-4, 1].

    3. Degree Estimation: The problem is now to separate two adjacent intervals [0, epsilon] and [epsilon, 1] where epsilon = n^-4.
       From approximation theory, the degree 'd' required for a polynomial to create a sharp transition in a region of relative size epsilon is d = Omega(1/sqrt(epsilon)).
       So, the degree of the polynomial in u, which is deg(q), is:
       deg(q) = Omega(1 / sqrt(n^-4)) = Omega(n^2).

    4. Final Result:
       We have d_n = deg(q) / 2.
       So, d_n = Omega(n^2).
       The asymptotic growth rate is d_n = Theta(n^alpha), so alpha must be 2.
    """
    alpha = 2
    print(f"The problem is transformed by the change of variable x = y^2.")
    print(f"This converts the problem for p_n(x) on [1, n^2] and [n^2+1, n^10] into a problem for an even polynomial q(y) on the adjacent intervals [1, n] and [n, n^5].")
    print(f"The degree of p_n(x) is half the degree of q(y).")
    print(f"Mapping the combined interval [1, n^5] to [0, 1], the two intervals for q(y) become approximately [0, n^-4] and [n^-4, 1].")
    print(f"The degree of a polynomial required to distinguish two adjacent intervals of relative size epsilon is proportional to 1/sqrt(epsilon).")
    print(f"Here, epsilon = n^-4.")
    print(f"So, the degree of q(y) is d_q = Theta(1/sqrt(n^-4)) = Theta(n^2).")
    print(f"The degree of p_n(x) is d_n = d_q / 2 = Theta(n^2).")
    print(f"Comparing with d_n = Theta(n^alpha), we find alpha = 2.")
    print(f"The value of alpha is {alpha}.")

solve()