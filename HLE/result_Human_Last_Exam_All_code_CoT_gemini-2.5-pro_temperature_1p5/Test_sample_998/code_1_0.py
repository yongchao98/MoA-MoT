import math

def sum_arithmetic_series(first, last, num_terms):
    """Calculates the sum of an arithmetic series."""
    return num_terms * (first + last) // 2

def main():
    """
    Solves the problem by calculating the components of the sum S,
    finding the prime factorization of S, and then computing the number of its divisors.
    """
    # Part 1: The constant factor from τ(2^8)
    part1 = 8 + 1

    # Part 2: The sum related to the prime 29 (exponent 59)
    # Sum_{y=0 to 59} (60-y) = 60+59+...+1
    part2 = sum_arithmetic_series(60, 1, 60)

    # Part 3: The double summation related to primes 59 (exp 79) and 79 (exp 29)
    # Sums for z from 0 to 79 (term is 80-z)
    sum_z_even = sum_arithmetic_series(80, 2, 40) # z=0,2,...,78 (40 terms)
    sum_z_odd = sum_arithmetic_series(79, 1, 40) # z=1,3,...,79 (40 terms)

    # Sums for w from 0 to 29 (term is 30-w)
    sum_w_even = sum_arithmetic_series(30, 2, 15) # w=0,2,...,28 (15 terms)
    sum_w_odd = sum_arithmetic_series(29, 1, 15) # w=1,3,...,29 (15 terms)

    part3 = sum_z_even * sum_w_even + sum_z_odd * sum_w_odd

    # The sum S is the product of these three parts.
    # S = part1 * part2 * part3

    # We need the prime factorization of S to find its number of divisors.
    # We find the exponents of the prime factors for each part.
    # part1 = 9 = 3^2
    # part2 = 1830 = 2 * 3 * 5 * 61
    # part3 = 753600 = 2^6 * 3 * 5^2 * 157
    
    # Exponents in the prime factorization of S = part1 * part2 * part3
    exponents = {
        2: 1 + 6,
        3: 2 + 1 + 1,
        5: 1 + 2,
        61: 1,
        157: 1
    }

    # Calculate the number of divisors (τ function)
    num_divisors = 1
    equation_parts = []
    
    # Sort keys for consistent output order
    sorted_primes = sorted(exponents.keys())

    for prime in sorted_primes:
        exp = exponents[prime]
        term = exp + 1
        num_divisors *= term
        equation_parts.append(f"({exp}+1)")

    # Print the equation for the number of divisors
    equation = " * ".join(equation_parts)
    print(f"The prime factorization of S leads to the exponents: {list(exponents.values())}")
    print(f"The number of divisors is the product of (exponent + 1) for each prime factor.")
    final_equation_values = " * ".join(str(exponents[p]+1) for p in sorted_primes)
    print(f"Number of divisors = {equation} = {final_equation_values} = {num_divisors}")
    print(f"The final answer is {num_divisors}.")

if __name__ == "__main__":
    main()
