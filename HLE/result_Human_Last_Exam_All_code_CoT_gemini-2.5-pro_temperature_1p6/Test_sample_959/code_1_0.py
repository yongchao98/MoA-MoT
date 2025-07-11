import math

def phi(n):
    """
    Computes Euler's totient function.
    It counts the number of positive integers up to n that are relatively prime to n.
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError("Input must be a positive integer.")
    
    count = 0
    for i in range(1, n + 1):
        if math.gcd(i, n) == 1:
            count += 1
    return count

def main():
    """
    Calculates the required sum based on the reasoning provided.
    """
    # As established in the explanation, the problem simplifies to finding phi(31).
    n = 31
    
    # Calculate phi(31)
    # For a prime number p, phi(p) = p - 1.
    result = n - 1
    
    # The collection of central extensions E contains a single element E = C_31.
    # The sum is the order of the outer automorphism group of C_31.
    # sum = o(C_31) = |Out(C_31)|
    # Since C_31 is abelian, |Out(C_31)| = |Aut(C_31)| = phi(31).
    final_sum = result

    # Output the final equation and the result
    print("The final sum is determined by the equation:")
    print(f"Sum = \u03C6({n}) = {n} - 1 = {final_sum}")

if __name__ == "__main__":
    main()
