import math

def jacobi_neg_one(k):
    """
    Calculates the Jacobi symbol (-1/k) for a positive odd integer k.
    """
    if k <= 0 or k % 2 == 0:
        raise ValueError("k must be a positive odd integer for Jacobi symbol (-1/k)")
    if k % 4 == 1:
        return 1
    else: # k % 4 == 3
        return -1

def calculate_zeta_n(n, N1, N2):
    """
    Calculates the higher central charge zeta_n for the U(1)_{2N1} x U(1)_{-2N2} theory.
    This formula is valid for odd n.
    """
    if n <= 0 or n % 2 == 0:
        print("This formula for zeta_n is valid only for positive odd integers n.")
        return

    # Calculate g1 = gcd(n, N1) and g2 = gcd(n, N2)
    g1 = math.gcd(n, N1)
    g2 = math.gcd(n, N2)

    # Calculate k1 = n/g1 and k2 = n/g2
    k1 = n // g1
    k2 = n // g2

    # Calculate the two Jacobi symbols
    jacobi1 = jacobi_neg_one(k1)
    jacobi2 = jacobi_neg_one(k2)

    # Calculate zeta_n
    zeta_n = jacobi1 * jacobi2

    # Print the results
    print(f"For n = {n}, N1 = {N1}, N2 = {N2}:")
    print(f"g1 = gcd(n, N1) = gcd({n}, {N1}) = {g1}")
    print(f"g2 = gcd(n, N2) = gcd({n}, {N2}) = {g2}")
    print(f"k1 = n / g1 = {n} / {g1} = {k1}")
    print(f"k2 = n / g2 = {n} / {g2} = {k2}")
    print(f"The formula for zeta_n is ((-1/k1)) * ((-1/k2))")
    print(f"zeta_n = ((-1/{k1})) * ((-1/{k2})) = ({jacobi1}) * ({jacobi2}) = {zeta_n}")


# Example usage:
n_val = 3
N1_val = 15
N2_val = 5
calculate_zeta_n(n_val, N1_val, N2_val)

print("\n" + "="*20 + "\n")

# Another example
n_val = 5
N1_val = 3
N2_val = 7
calculate_zeta_n(n_val, N1_val, N2_val)