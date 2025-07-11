import math

def solve_crawley_nordstrom():
    """
    Calculates the Crawley-Nordström invariant for a given Calabi-Yau Link.
    """
    # Step 1: Define weights and degree.
    weights = [22, 29, 49, 50, 75]
    
    # As explained in the text, we assume the degree d=225.
    d = 225

    # Step 2: Calculate the components of the formula.
    d_squared = d**2
    weights_squared = [w**2 for w in weights]
    sum_weights_squared = sum(weights_squared)

    # Step 3: Calculate the invariant.
    numerator = d_squared - sum_weights_squared
    denominator = d
    c_invariant = numerator / denominator

    # Step 4: Print the detailed equation and the result, as requested.
    print(f"The Crawley-Nordström invariant 'c' is calculated using the formula:")
    print(f"c = (d^2 - Σ(wi^2)) / d")
    print(f"\nGiven weights: {weights}")
    print(f"Assumed degree d: {d}")

    print(f"\nCalculation steps:")
    print(f"d^2 = {d}^2 = {d_squared}")

    w_sq_str = ' + '.join([str(w) + '^2' for w in weights])
    w_sq_val_str = ' + '.join(map(str, weights_squared))
    print(f"Σ(wi^2) = {w_sq_str}")
    print(f"         = {w_sq_val_str}")
    print(f"         = {sum_weights_squared}")

    print("\nSubstituting the values into the formula:")
    print(f"c = ({d_squared} - {sum_weights_squared}) / {d}")
    print(f"c = {numerator} / {denominator}")
    
    # Check for common factors to present a simplified fraction
    common_divisor = math.gcd(numerator, denominator)
    if common_divisor > 1:
        simplified_num = numerator // common_divisor
        simplified_den = denominator // common_divisor
        print(f"c = {simplified_num} / {simplified_den} (in simplest form)")
    
    print("\nFinal Result:")
    print(f"The Crawley-Nordström invariant is {c_invariant}")


solve_crawley_nordstrom()