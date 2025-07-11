import sys

def calculate_chi_star(chi):
    """
    Calculates chi* based on the relationship chi* = -chi / (1 + chi).

    Args:
        chi (float): The magnetic susceptibility of the material.

    Returns:
        float or None: The calculated value of chi* or None if the calculation is undefined.
    """
    # Check for the undefined case where the denominator is zero
    if 1 + chi == 0:
        print("Error: The value of chi cannot be -1, as it leads to a division by zero.", file=sys.stderr)
        return None

    # Calculate chi* using the derived formula
    chi_star = -chi / (1 + chi)
    return chi_star

def main():
    """
    Main function to demonstrate the calculation of chi*.
    """
    # Example value for chi
    chi = 3.0
    
    print(f"Given the magnetic susceptibility chi = {chi}")
    print("The relationship between chi and chi* is: chi* = -chi / (1 + chi)")
    
    # Calculate chi_star
    chi_star = calculate_chi_star(chi)
    
    if chi_star is not None:
        # Print the final equation with the numbers substituted
        print("\nSubstituting the value of chi into the equation:")
        print(f"chi* = -{chi} / (1 + {chi})")
        print(f"chi* = {chi_star:.4f}")

if __name__ == "__main__":
    main()