import math

def solve():
    """
    Solves the problem based on the analysis that for large n, l(a) simplifies to -2/a.
    The contradictions in the problem statement regarding the matrix Sigma suggest that
    the term involving Sigma must vanish.
    """

    # The first 10 prime numbers
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]

    total_sum = 0
    l_values = []
    for p in primes:
        # Based on the derivation, l(a) simplifies to -2/a for large n.
        l_val = -2.0 / p
        l_values.append(l_val)
        total_sum += l_val

    final_answer = math.floor(total_sum)

    # Building the equation string for the final output, as requested.
    equation_parts = []
    for i in range(len(primes)):
        p = primes[i]
        val = l_values[i]
        # Format to show a few decimal places for clarity
        equation_parts.append(f"{val:.8f}")

    sum_str = " + ".join(equation_parts)
    
    print(f"The calculation is the floor of the sum of l(a_i) for the first 10 primes a_i.")
    print(f"l(a) is approximated as -2/a.")
    print(f"The sum is floor(l(2) + l(3) + ... + l(29))")
    print(f"Final equation: floor({sum_str}) = {final_answer}")


solve()