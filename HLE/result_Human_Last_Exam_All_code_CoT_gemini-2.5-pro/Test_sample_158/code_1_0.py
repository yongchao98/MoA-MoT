def solve_magnetic_cloaking():
    """
    This function explains and prints the solution to the magnetic cloaking problem.

    The solution is derived based on the physical assumption that the interior
    region (rho < R1) is a perfect magnetic conductor (PMC), where permeability mu -> infinity.
    This assumption is necessary to find a unique, non-trivial solution as requested
    by the problem statement.

    A PMC core forces the magnetic field H inside it to be zero. The specific
    permeability of the shell then works to cancel the distortion the PMC core
    would otherwise cause in the exterior region.
    """

    # Explanation of the results
    print("To find a non-trivial solution, we assume the interior region (rho < R1) is a perfect magnetic conductor.")
    print("This leads to a specific permeability for the shell that cloaks the core, leaving the exterior field undisturbed.")
    print("The consequence of this assumption is that the magnetic field H inside the core is zero.")
    print("-" * 50)

    # Result for the permeability of the shell
    print("The required permeability 'mu' of the magnetic material in the shell is given by the formula:")
    # The formula is printed as a string to ensure exact formatting, including all numbers.
    # The equation is: mu = mu_0 * (R2^2 - R1^2) / (R2^2 + R1^2)
    print("mu = mu_0 * (R2**2 - R1**2) / (R2**2 + R1**2)")
    print("\nwhere:")
    print("  mu_0 is the permeability of free space.")
    print("  R1 is the interior radius of the shell.")
    print("  R2 is the exterior radius of the shell.")
    print("-" * 50)

    # Result for the magnetic field in the interior region
    print("The magnetic field 'H_int' in the interior region (rho < R1) is:")
    # The equation is: H_int = 0
    print("H_int = 0")
    print("-" * 50)

# Execute the function to print the solution
solve_magnetic_cloaking()
