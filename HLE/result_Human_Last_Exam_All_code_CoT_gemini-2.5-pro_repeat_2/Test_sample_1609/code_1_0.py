import sys

def get_max_prime_implicants(n):
    """
    Provides the maximal number of prime implicants for a Boolean function of n variables.

    This function uses pre-computed, known values from the OEIS sequence A001035,
    as direct computation is extremely complex and computationally expensive.
    """
    # a(n): Maximal number of prime implicants of a Boolean function of n variables.
    # This is the sequence A001035 in the On-Line Encyclopedia of Integer Sequences (OEIS).
    known_values = {
        0: 1,
        1: 2,
        2: 6,
        3: 20,
        4: 78,
        5: 352,
        6: 1792,
    }

    if n in known_values:
        return known_values[n]
    else:
        # For n > 6, the values are known but get very large.
        # This function is limited to its pre-computed table.
        return None

def main():
    """
    Main function to solve the problem for n=4.
    """
    n = 4
    result = get_max_prime_implicants(n)

    if result is not None:
        # The final "equation" is a(n) = result.
        # We output the numbers n and result as requested.
        print(f"The maximal number of prime implicants for a Boolean function of n = {n} variables is:")
        print(f"a({n}) = {result}")
    else:
        print(f"The value for a({n}) is not in our pre-computed list.", file=sys.stderr)

if __name__ == "__main__":
    main()