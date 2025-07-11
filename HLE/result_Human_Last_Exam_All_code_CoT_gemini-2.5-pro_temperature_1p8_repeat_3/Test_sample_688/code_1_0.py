import math

def calculate_cn(n_max):
    """
    Calculates and prints the prefactor c_n for the virial expansion
    from n=2 to n_max.

    The formula for c_n is:
    c_n = - (1/n) * (-1)^(((n-2)*(n-3))/2)
    """
    print(f"Calculating the prefactor c_n for n from 2 to {n_max}")
    print("-" * 30)

    for n in range(2, n_max + 1):
        # Calculate the exponent for the sign term
        exp_term1 = n - 2
        exp_term2 = n - 3
        
        # Using integer division
        exponent = (exp_term1 * exp_term2) // 2

        # Calculate the final numerator and denominator
        sign = (-1)**exponent
        numerator = -1 * sign
        denominator = n
        
        # Simplify fraction for printing
        common_divisor = math.gcd(numerator, denominator)
        num_simple = numerator // common_divisor
        den_simple = denominator // common_divisor

        # Print the detailed calculation
        print(f"For n = {n}:")
        print(f"  c_{n} = - (1/{n}) * (-1)^((({n}-2)*({n}-3))/2)")
        print(f"  c_{n} = - (1/{n}) * (-1)^(({exp_term1}*{exp_term2})/2)")
        print(f"  c_{n} = - (1/{n}) * (-1)^{exponent}")
        print(f"  c_{n} = - (1/{n}) * {sign}")
        print(f"  c_{n} = {numerator}/{denominator}")

        if den_simple == 1:
            print(f"  Final value: c_{n} = {num_simple}")
        else:
            print(f"  Final value: c_{n} = {num_simple}/{den_simple} = {num_simple/den_simple:.4f}")
        
        print("-" * 30)

# Calculate c_n for n from 2 to 10
calculate_cn(10)
