import math

def calculate_cn(n):
    """
    Calculates and prints the derivation for the prefactor c_n for the
    fully f-connected Ree-Hoover diagram.

    Args:
        n (int): The number of particles in the diagram (must be >= 2).
    """
    if not isinstance(n, int) or n < 2:
        print(f"Error: n must be an integer >= 2. Got {n}.")
        return

    print("-" * 40)
    print(f"Calculating the prefactor c_n for n = {n}:")

    # General formula for the contribution of a single Ree-Hoover diagram
    print("\nThe contribution of a single Ree-Hoover diagram 'g' to the virial coefficient B_n is:")
    print("  Contribution = - (w_g / n) * Lambda_g")
    print("where w_g is the diagram weight and Lambda_g is the cluster integral.")
    print("The prefactor c_n is therefore given by c_n = -w_g / n.")

    if n == 2:
        weight = 1
        c_n_num = -1
        c_n_den = 2
        
        print(f"\nFor the fully f-connected diagram with n = {n}:")
        print(f"The Ree-Hoover weight 'w' is {weight}.")
        print(f"The prefactor c_{n} is therefore:")
        print(f"  c_{n} = - (w / n) = - ({weight}) / {n}")
        print(f"  c_{n} = {c_n_num}/{c_n_den} = {c_n_num/c_n_den}")
        
    else:  # n >= 3
        weight_val = -math.factorial(n - 2)
        c_n_num = math.factorial(n - 2)
        c_n_den = n
        
        print(f"\nFor the fully f-connected diagram with n = {n}:")
        print(f"The Ree-Hoover weight 'w' is -(n-2)! = -({n}-2)! = -{math.factorial(n-2)} = {weight_val}")
        print(f"The prefactor c_{n} is therefore:")
        print(f"  c_{n} = - (w / n) = - ({weight_val}) / {n} = {c_n_num} / {n}")

        # Simplify the fraction if possible by finding the greatest common divisor
        common_divisor = math.gcd(c_n_num, c_n_den)
        if common_divisor > 1:
            print(f"  As a simplified fraction: c_{n} = {c_n_num//common_divisor}/{c_n_den//common_divisor}")
        
        print(f"  As a decimal: c_{n} = {c_n_num/c_n_den}")

if __name__ == "__main__":
    # Since the user did not specify a value for n, we will demonstrate the
    # calculation for a few small values. You can change the values in the
    # list below to compute c_n for any integer n >= 2.
    n_values = [2, 3, 4, 5]
    
    for n in n_values:
        calculate_cn(n)