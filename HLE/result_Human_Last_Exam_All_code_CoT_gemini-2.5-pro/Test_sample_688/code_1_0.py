import math

def calculate_virial_prefactor(n):
    """
    Calculates the system-independent prefactor c_n for the fully f-connected
    diagram's contribution to the n-th virial coefficient.
    
    Args:
        n (int): The order of the virial coefficient (must be an integer >= 2).
    """
    if not isinstance(n, int) or n < 2:
        print("Error: The order of the virial coefficient, n, must be an integer greater than or equal to 2.")
        return

    # The formula for the prefactor is c_n = -(n-1) / n!
    # This is derived from the standard Mayer cluster expansion for virial coefficients.
    
    # Calculate the numerator of the fraction.
    numerator = n - 1
    
    # Calculate the denominator of the fraction (n factorial).
    denominator = math.factorial(n)
    
    # Calculate the final value of the prefactor.
    c_n_value = -numerator / denominator
    
    # The final output should include each number in the equation, as requested.
    # We construct a string to show the full calculation.
    print(f"The formula for the prefactor is: c_n = -(n-1) / n!")
    print(f"For n = {n}, the calculation is as follows:")
    print(f"c_{n} = -({n} - 1) / {n}! = -{numerator} / {denominator} = {c_n_value}")

# --- Main execution ---
# You can change the value of n below to calculate the prefactor for a different
# virial coefficient.
n_value = 4
calculate_virial_prefactor(n_value)