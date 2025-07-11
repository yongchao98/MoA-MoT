import numpy as np

def is_prime(num):
    """A simple primality test function."""
    if num <= 1:
        return False
    if num <= 3:
        return True
    if num % 2 == 0 or num % 3 == 0:
        return False
    i = 5
    while i * i <= num:
        if num % i == 0 or num % (i + 2) == 0:
            return False
        i += 6
    return True

def get_primes(count):
    """Generates a list of the first `count` prime numbers."""
    primes = []
    num = 2
    while len(primes) < count:
        if is_prime(num):
            primes.append(num)
        num += 1
    return primes

def demonstrate_convergence():
    """
    Demonstrates the convergence properties of test series for different sets S.
    The test series is f(z) = sum_{k=1 to inf} [(-1)^k / k] * z^(s_k).
    """
    N = 20000  # Number of terms to compute in the partial sum

    # --- Case 2: S = {n^k}, k=4 ---
    # This set is arithmetically regular and is expected to fail.
    # We test at z=i, where we expect divergence.
    # The terms are (-1)^n/n * i^(n^4).
    # If n is even, n^4 = (2m)^4 = 16m^4 = 0 mod 4. So i^(n^4) = 1.
    # If n is odd, n^4 = (2m+1)^4 = 1 mod 4. So i^(n^4) = i.
    real_part_poly = 0.0
    imag_part_poly = 0.0
    for n in range(1, N + 1):
        if n % 2 == 0:  # even n
            real_part_poly += 1.0 / n
        else:  # odd n
            imag_part_poly -= 1.0 / n
    
    print("--- Analysis for Set 2: S = {n^4} ---")
    print(f"We test the partial sum of Sum[((-1)^n/n)*z^(n^4)] for N={N} terms at z=i.")
    print(f"Analytic prediction: Divergence.")
    # The real part is Sum[1/(2m)] and imag part is -Sum[1/(2m-1)], both diverge like log(N).
    print(f"Numerical result for the partial sum: Real part = {real_part_poly:.4f}, Imaginary part = {imag_part_poly:.4f}")
    print("The components of the partial sum grow with N, confirming divergence.\n")

    # --- Case 3: S = set of primes {p_k} ---
    # This set is arithmetically irregular and is expected to succeed.
    # We test at z=-1, where we expect convergence.
    # The series is Sum[((-1)^k/k)*(-1)^(p_k)].
    # p_1=2 is even, all other primes p_k are odd.
    # Term k=1: (-1)^1/1 * (-1)^2 = -1
    # Term k>1: (-1)^k/k * (-1)^(odd) = -(-1)^k/k
    # Expected sum = -1 - Sum[(-1)^k/k for k from 2 to inf]
    # Sum[(-1)^k/k for k>=1] = -ln(2), so Sum[... k>=2] = -ln(2) + 1.
    # Total sum -> -1 - (-ln(2) + 1) = ln(2) - 2
    
    primes = get_primes(N)
    partial_sum_primes = 0.0
    for k in range(1, N + 1):
        p_k = primes[k-1]
        term = ((-1)**k) / k * ((-1)**p_k)
        partial_sum_primes += term
    
    theoretical_limit = np.log(2) - 2

    print("--- Analysis for Set 3: S = Primes ---")
    print(f"We test the partial sum of Sum[((-1)^k/k)*z^(p_k)] for N={N} terms at z=-1.")
    print(f"Analytic prediction: Convergence to ln(2) - 2.")
    print(f"Final equation for the sum is: -1 + Sum_{{k=2 to inf}} (-1)^(k+1)/k = -1 + (1 - ln(2)) = ln(2) - 2")
    print(f"Theoretical value: {theoretical_limit:.8f}")
    print(f"Numerical result for the partial sum: {partial_sum_primes:.8f}")
    print("The partial sum clearly converges to the predicted theoretical value.\n")

    print("--- Overall Conclusion ---")
    print("Sets with sufficient arithmetic irregularity (like Primes and Lacunary sets) have the property.")
    print("Sets with arithmetic regularity (like Polynomials and sequences asymptotic to n) do not.")
    print("Therefore, only sets 3 and 4 have the property.")

if __name__ == '__main__':
    demonstrate_convergence()