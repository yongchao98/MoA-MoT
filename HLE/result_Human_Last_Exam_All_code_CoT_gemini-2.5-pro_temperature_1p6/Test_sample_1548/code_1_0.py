import math

def demonstrate_infinite_transient_set():
    """
    This function demonstrates that not every transient set for the specified
    Markov chain must be finite. It does so by constructing an infinite set
    and showing that it meets the criterion for being transient.
    """
    print("The statement 'every transient set must necessarily be finite' is false.")
    print("A set A is transient for the h-transformed simple random walk in Z^2 if and only if the following sum converges:")
    print("  Sum_{x in A} (log ||x||)^(-2) < infinity")
    print("\nWe can construct an infinite set for which this sum is finite, proving the statement is false.")
    print("Consider the infinite set A = {x_k = (r_k, 0) | k = 1, 2, 3, ...}, where the norm r_k = 2^(2^k).")
    print("Let's check if this set is transient by computing the sum.\n")

    # The term for x_k in the series is 1 / (log(r_k))^2.
    # log(r_k) = log(2^(2^k)) = 2^k * log(2).
    # So the term is 1 / (2^k * log(2))^2 = (1 / (log(2))^2) * (1 / (4^k)).
    # The sum is a geometric series multiplied by a constant.
    print("The k-th term of the series is: 1 / (log(2^(2^k)))^2")
    print("\n--- Numerical Demonstration of Convergence ---")
    print(f"{'k':<5}{'Norm ||x_k||':<15}{'Term Value':<25}{'Partial Sum':<25}")
    print("-" * 75)
    
    partial_sum = 0
    max_k = 10
    for k in range(1, max_k + 1):
        # We represent the norm symbolically to avoid huge numbers
        norm_str = f"2^(2^{k})"
        # The term can be calculated without evaluating the large norm
        term = 1 / ((2**k * math.log(2))**2)
        partial_sum += term
        print(f"{k:<5}{norm_str:<15}{term:<25.15e}{partial_sum:<25.15f}")
    
    print("\nThe partial sum clearly converges to a finite value.")
    
    print("\n--- Analytical Calculation of the Final Equation ---")
    # The infinite series is Sum_{k=1 to inf} (1 / log(2)^2) * (1/4)^k
    # The sum of the geometric series Sum_{k=1 to inf} (1/4)^k is (1/4) / (1 - 1/4) = 1/3.
    log2 = math.log(2)
    series_sum = (1/3)
    constant_factor = 1 / (log2**2)
    analytical_sum = constant_factor * series_sum

    print("The final equation for the sum is:")
    print(f"  Sum = [ Sum_{{k=1 to inf}} (1/4)^k ] / (log(2))^2")
    print("\nPlugging in the numbers for each part:")
    print(f"  Sum of the geometric series = {series_sum:.6f}")
    print(f"  log(2) = {log2:.6f}")
    print(f"  Final Sum = {series_sum:.6f} / ({log2:.6f})^2 = {analytical_sum:.6f}")

    print("\nSince the sum is finite, the chosen infinite set A is a transient set.")
    print("Therefore, it is not true that every transient set must be finite.")

demonstrate_infinite_transient_set()
<<<No>>>