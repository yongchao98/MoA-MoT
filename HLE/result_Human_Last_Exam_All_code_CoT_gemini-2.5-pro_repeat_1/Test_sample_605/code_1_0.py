import sys

def calculate_crawley_nordstrom_invariant():
    """
    Calculates the Crawley-Nordstrom invariant for a Calabi-Yau hypersurface.

    The Crawley-Nordstrom invariant is defined as c_N = d - sum(w_i), where 'd' is the
    degree of the defining quasi-homogeneous polynomial and 'w_i' are the weights of the
    ambient weighted projective space.

    For a Calabi-Yau hypersurface, the condition d = sum(w_i) must hold.
    This means the invariant is always 0 for a (well-defined) Calabi-Yau variety.
    """
    # The weights for the variables (z_1, z_2, z_3, z_4, z_5)
    weights = [22, 29, 49, 50, 75]

    # Calculate the sum of the weights
    sum_of_weights = sum(weights)

    # For a Calabi-Yau variety, the degree 'd' of the polynomial is equal to the sum of the weights.
    d = sum_of_weights

    # Calculate the Crawley-Nordstrom invariant
    invariant = d - sum_of_weights
    
    # Construct the string for the sum of weights part of the equation
    weights_sum_str = " + ".join(map(str, weights))

    print(f"The weights are: w = {weights}")
    print(f"The sum of the weights is: {weights_sum_str} = {sum_of_weights}")
    print(f"For a Calabi-Yau variety, the degree 'd' of the polynomial must equal the sum of the weights.")
    print(f"Therefore, the degree d = {d}.")
    print("\nThe Crawley-Nordström invariant is given by the formula: d - (sum of weights)")
    print(f"Final Equation: {d} - {sum_of_weights} = {invariant}")
    
    # Appending the final answer in the requested format to stdout
    sys.stdout.write(f"\n<<<{invariant}>>>")

calculate_crawley_nordstrom_invariant()