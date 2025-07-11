import math

def get_continued_fraction_coeffs(p, q):
    """
    Computes the coefficients of the simple continued fraction of p/q.
    Uses the Euclidean algorithm.
    """
    if q == 0:
        raise ValueError("Denominator cannot be zero.")
    
    # Ensure q is positive for consistent representation
    if q < 0:
        p, q = -p, -q
        
    coeffs = []
    while q != 0:
        a = p // q
        coeffs.append(a)
        p, q = q, p % q
    return coeffs

def main():
    """
    Computes and prints the continued fraction for the generalized Markov number m_{4/7}.
    """
    p, q = 4, 7

    print(f"Goal: Compute the continued fraction for the generalized Markov number m_{p}/{q}.")
    print("-" * 50)

    # Step 1: Compute the continued fraction of p/q
    cf_pq = get_continued_fraction_coeffs(p, q)
    cf_pq_str = f"[{cf_pq[0]}; " + ", ".join(map(str, cf_pq[1:])) + "]"
    print(f"Step 1: The continued fraction of {p}/{q} is {cf_pq_str}.")
    print("-" * 50)

    # Step 2: Construct the periodic part for m_{p/q}
    # The construction formula assumes 0 < p/q < 1, which means a_0 = 0.
    # This is true for 4/7.
    a0 = cf_pq[0]
    a = cf_pq[1:] # This is the sequence [a_1, a_2, ..., a_n]
    n = len(a)

    if n == 0: # Integer case, not applicable here
        periodic_part = []
    elif n == 1: # e.g., for 1/3 = [0;3] -> a=[3]
        periodic_part = [a[0] - 1, 1]
    else: # General case for n > 1
        first_part = a[:-1]
        middle_part = [a[-1] - 1, 1]
        last_part = list(reversed(first_part))
        periodic_part = first_part + middle_part + last_part

    print("Step 2: Construct the periodic part using the formula:")
    print("m_{p/q} = [a_0; overline(a_1, ..., a_{n-1}, a_n-1, 1, a_{n-1}, ..., a_1)]")
    print(f"For {p}/{q}, the coefficients [a_1, ..., a_n] are {a}.")
    print(f"The first part (a_1 to a_{n-1}) is {a[:-1]}.")
    print(f"The middle part (a_n-1, 1) is [{a[-1]-1}, 1].")
    print(f"The last part (reversed first part) is {list(reversed(a[:-1]))}.")
    print("-" * 50)
    
    # Step 3: Print the final result
    print("Step 3: The resulting continued fraction is:")
    
    # This loop prints each number as requested by the prompt.
    print(f"m_{p}/{q} = [{a0}; ", end="")
    
    # Use "overline" for terminal display of periodic part
    periodic_part_str = ", ".join(map(str, periodic_part))
    print(f"overline({periodic_part_str})]")

if __name__ == "__main__":
    main()