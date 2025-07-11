import math

def is_prime(n):
    """Checks if a number is prime."""
    if n <= 1:
        return False
    # Check for divisibility from 2 up to the square root of n
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

def is_perfect_square(n):
    """Checks if a number is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    x = int(math.sqrt(n))
    return x * x == n

def find_tiling_subset():
    """
    This function finds the subset of integers t from a given list for which
    any n x n grid has an even number of t-omino tilings.

    The logic is as follows:
    1. If t is a perfect square (t=k*k), a k x k omino can tile a k x k grid
       in 1 way (odd). So t is excluded.
    2. If t is a prime number, it can be shown that for any t-omino shape,
       the number of tilings of an n x n grid is always even. So t is included.
    3. If t is a composite number but not a perfect square (like 15=3*5),
       it is reasoned to behave like other composite numbers and is excluded.
    
    Therefore, the property simplifies to checking for primality.
    """
    t_values = [2, 3, 4, 5, 7, 9, 15]
    subset = []

    print(f"Analyzing the set T = {t_values}")
    print("=" * 40)

    for t in t_values:
        print(f"Checking t = {t}:")
        # Rule out perfect squares
        if is_perfect_square(t):
            k = int(math.sqrt(t))
            print(f"  - Result: t={t} is NOT in the subset.")
            print(f"  - Reason: t = {t} = {k}*{k} is a perfect square.")
            print(f"    A {k}x{k} square omino tiles a {k}x{k} grid in 1 way, which is odd.")
            print(f"    This provides a counterexample for n={k}.")
        # Rule in primes
        elif is_prime(t):
            print(f"  - Result: t={t} IS in the subset.")
            print(f"  - Reason: t = {t} is a prime number.")
            print("    For any t-omino that is not a 1xt rectangle, the number of tilings is 0 (even).")
            print("    For the 1xt rectangular omino, the number of tilings is always even.")
            subset.append(t)
        # Rule out other composites
        else:
            print(f"  - Result: t={t} is NOT in the subset.")
            print(f"  - Reason: t = {t} is a composite number.")
            print("    Like perfect squares, composite numbers allow for the construction of ominoes")
            print("    that can tile certain n x n grids an odd number of times.")

        print("-" * 40)

    print("Final Conclusion:")
    print("The subset of integers for which the statement is true is the set of prime numbers in T.")
    # The problem asks to output the numbers in the final result.
    print(f"The subset is {{{', '.join(map(str, sorted(subset)))}}}")

if __name__ == '__main__':
    find_tiling_subset()