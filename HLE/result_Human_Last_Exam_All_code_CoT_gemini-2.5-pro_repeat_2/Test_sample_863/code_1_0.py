import sys

def solve_for_chi_star(a, b, chi):
    """
    Calculates the susceptibility chi* based on the given equation.

    Args:
        a (float): Half-dimension of the prism along the x-direction.
        b (float): Half-dimension of the prism along the y-direction.
        chi (float): Magnetic susceptibility of the material.

    Returns:
        float: The calculated susceptibility chi*.
    """
    if a <= 0 or b <= 0:
        print("Error: Dimensions 'a' and 'b' must be positive.")
        sys.exit(1)

    # Step 1: Calculate the geometric demagnetizing factors Nx and Ny
    Nx = b / (a + b)
    Ny = a / (a + b)

    print(f"Geometric factors: Nx = {Nx:.4f}, Ny = {Ny:.4f}")
    print(f"Check: Nx + Ny = {Nx + Ny:.4f}\n")

    # Step 2: Calculate the first term, Nm(a/b, chi)
    nm_1 = Nx / (1 + Nx * chi)

    print(f"Given a = {a}, b = {b}, and chi = {chi}:")
    print(f"The first term is Nm(a/b, chi) = {nm_1:.4f}\n")

    # Step 3: From the main equation, find the required value for the second term
    # Nm(a/b, chi) + Nm(b/a, chi*) = 1  =>  Nm(b/a, chi*) = 1 - Nm(a/b, chi)
    nm_2_target = 1 - nm_1
    
    print(f"For the equation to hold, the second term Nm(b/a, chi*) must be {nm_2_target:.4f}\n")

    # Step 4: Solve for chi* using the definition of the second term
    # nm_2 = Ny / (1 + Ny * chi*)
    # 1 + Ny * chi* = Ny / nm_2
    # Ny * chi* = (Ny / nm_2) - 1
    # chi* = (1 / nm_2) - (1 / Ny)
    if nm_2_target == 0:
        print("Error: Cannot solve for chi* because Nm(b/a, chi*) is zero, which implies infinite chi*.")
        return float('inf')
        
    chi_star = (1 / nm_2_target) - (1 / Ny)

    # Step 5: Verification
    # Recalculate the second term using the found chi*
    nm_2_calculated = Ny / (1 + Ny * chi_star)
    
    # The final equation with the calculated numbers
    final_sum = nm_1 + nm_2_calculated

    print("--- Verification ---")
    print(f"Calculated chi* = {chi_star:.4f}")
    print("Plugging chi* back into the equation for the second term:")
    print(f"Nm(b/a, chi*) = Ny / (1 + Ny * chi*) = {nm_2_calculated:.4f}")
    print("\nFinal equation with calculated values:")
    print(f"Nm(a/b, chi) + Nm(b/a, chi*) = {nm_1:.4f} + {nm_2_calculated:.4f} = {final_sum:.4f}")

    return chi_star

if __name__ == '__main__':
    # --- User-defined parameters ---
    # You can change these values to explore different scenarios.
    
    # a and b are the semi-axes of the elliptical cross-section
    a = 2.0
    b = 1.0
    
    # chi is the intrinsic magnetic susceptibility
    chi = 5.0
    # ------------------------------------
    
    chi_star_result = solve_for_chi_star(a, b, chi)

    print(f"\nFor a={a}, b={b}, and chi={chi}, the required susceptibility is chi* = {chi_star_result:.4f}")
    # The final answer is returned in the specified format below for the example values.
    # e.g., <<< -0.3571 >>>
    final_answer = f"<<< {chi_star_result:.4f} >>>"