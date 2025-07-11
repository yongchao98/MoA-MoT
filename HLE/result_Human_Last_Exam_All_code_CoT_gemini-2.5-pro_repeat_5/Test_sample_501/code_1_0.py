import math

def display_force_law():
    """
    This function prints the derived force law for a thermally isolated polymer chain.
    """
    
    # The derived force law is of the form F(x) = -k * x
    # F(x) = - (2 * E(0) * x) / (n^2 * l^2)
    
    # Define the numerical constants and powers in the equation
    coefficient = 2
    power_of_n = 2
    power_of_l = 2
    
    print("The force law between the ends of a thermally isolated polymer is a restoring force, F(x) = -k*x.")
    print("This is valid for small end-to-end separation, x.")
    print("\nThe detailed equation for the force, F(x), is:")
    
    # The use of f-strings here helps format the output clearly, showing each component.
    print(f"F(x) = - ( {coefficient} * E(0) * x ) / ( n^{power_of_n} * l^{power_of_l} )")
    
    print("\nWhere the terms are:")
    print("--------------------")
    print(f"Numerical Coefficient: {coefficient}")
    print(f"E(0): The initial kinetic energy of the polymer at zero extension (x=0).")
    print(f"x:    The separation of the polymer ends.")
    print(f"n:    The number of links in the polymer chain (raised to the power {power_of_n}).")
    print(f"l:    The length of each link (raised to the power {power_of_l}).")

if __name__ == '__main__':
    display_force_law()
