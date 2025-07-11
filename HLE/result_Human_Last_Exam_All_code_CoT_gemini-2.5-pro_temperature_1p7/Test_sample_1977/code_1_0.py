import math

def get_norm(n):
    """
    Calculates the 1-norm of the correlation matrix T for the state J_n for a given even integer n.
    """
    if n % 2 != 0:
        raise ValueError("n must be an even integer.")

    # n_plus_1 is the number of pairs, and is odd
    n_plus_1 = n + 1
    
    total_sum_val = 0
    
    # Iterate over all possible counts of k0, k1, k2, k3
    # k0 is the count of sigma_0 (Identity)
    # k1 is the count of sigma_1 (Pauli-X)
    # k2 is the count of sigma_2 (Pauli-Y)
    # k3 is the count of sigma_3 (Pauli-Z)
    for k0 in range(n_plus_1 + 1):
        for k1 in range(n_plus_1 - k0 + 1):
            for k2 in range(n_plus_1 - k0 - k1 + 1):
                k3 = n_plus_1 - k0 - k1 - k2
                if k3 < 0:
                    continue

                # The correlation matrix T is indexed by non-trivial local operators,
                # so the case where all operators are Identity (k0 = n+1) is excluded.
                if k0 == n_plus_1:
                    continue
                
                # Number of ways to arrange these operators
                num_combinations = math.comb(n_plus_1, k0) * math.comb(n_plus_1 - k0, k1) * math.comb(n_plus_1 - k0 - k1, k2)

                # The coefficient depends on the parities of the counts
                # c = (1, 1, -1, 1), d = (3, -1, 1, -1)
                # product of c's = (-1)^k2
                # product of d's = 3^k0 * (-1)^k1 * (-1)^k3 = 3^k0 * (-1)^(k1+k3)
                
                # t_mu is proportional to (-1)^k2 + (1/3) * 3^k0 * (-1)^(k1+k3)
                # Since n+1 is odd, k1+k3 = n+1 - k0 - k2 is odd if k0+k2 is even, and even if k0+k2 is odd.
                # So (-1)^(k1+k3) = -(-1)^(k0+k2)
                
                # v = (-1)^k2 - 3^(k0-1) * (-1)^(k0+k2)
                # |v| = |1 - 3^(k0-1)*(-1)^k0|
                if k0 == 0:
                    abs_v = 2.0 / 3.0
                elif k0 % 2 == 1: # k0 is odd
                    abs_v = 1 + 3**(k0 - 1)
                else: # k0 is even and >= 2
                    abs_v = 3**(k0 - 1) - 1
                
                total_sum_val += num_combinations * abs_v

    norm_T1 = total_sum_val / (1 + 3**n)
    return round(norm_T1)

# Main part of the script
if __name__ == '__main__':
    # Calculate for n=0 and n=2 to demonstrate the pattern
    n_0 = 0
    norm_0 = get_norm(n_0)
    print(f"For n = {n_0}, the 1-norm is {norm_0}")
    print(f"The equation 3^n gives: 3^{n_0} = {3**n_0}")
    print("-" * 20)

    n_2 = 2
    norm_2 = get_norm(n_2)
    print(f"For n = {n_2}, the 1-norm is {norm_2}")
    print(f"The equation 3^n gives: 3^{n_2} = {3**n_2}")
    print("-" * 20)

    # The final equation
    print("Based on the pattern, the general equation for the 1-norm for even n is:")
    print("||T||_1 = 3^n")
