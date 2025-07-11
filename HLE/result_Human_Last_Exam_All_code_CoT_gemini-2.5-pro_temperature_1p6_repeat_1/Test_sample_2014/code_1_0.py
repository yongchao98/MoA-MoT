import sys

def solve_aerofoil_problem():
    """
    Calculates the lift ratio L1/L2 for two aerofoils in tandem formation
    in ground effect using the mirror image method.
    """
    # The problem can be solved by letting the chord c = 1, as the ratio
    # is dimensionless and c will cancel out.
    c = 1.0

    # Given parameters
    s = 0.5 * c  # Separation distance
    h = 0.5 * c  # Ride height

    print(f"Step 1: Define system parameters relative to chord c={c}.")
    print(f"Separation distance, s = {s:.1f}c")
    print(f"Ride height, h = {h:.1f}c\n")

    # The ratio of the lifts L1/L2 can be found with the formula:
    # L1 / L2 = (1 - C) / (1 + C)
    # where C is an interaction coefficient given by:
    # C = (c * 2 * h**2) / (s * (s**2 + 4 * h**2))

    print("Step 2: Calculate the interaction coefficient 'C'.")
    print("Formula: C = (c * 2 * h**2) / (s * (s**2 + 4 * h**2))")

    # Calculate numerator and denominator of C
    c_num_val = c * 2 * h**2
    c_den_val = s * (s**2 + 4 * h**2)
    C = c_num_val / c_den_val
    
    print(f"  Numerator = {c:.1f} * 2 * {h:.1f}**2 = {c_num_val:.3f}")
    print(f"  Denominator = {s:.1f} * ({s:.1f}**2 + 4 * {h:.1f}**2) = {c_den_val:.3f}")
    print(f"C = {c_num_val:.3f} / {c_den_val:.3f} = {C:.2f}\n")
    
    # Calculate the final lift ratio
    print("Step 3: Calculate the lift ratio L1/L2.")
    print("Formula: L1/L2 = (1 - C) / (1 + C)")

    # Numerator of the ratio
    ratio_num = 1 - C
    # Denominator of the ratio
    ratio_den = 1 + C
    # Final ratio
    lift_ratio = ratio_num / ratio_den
    
    print(f"  L1/L2 = (1 - {C:.2f}) / (1 + {C:.2f})")
    print(f"  L1/L2 = {ratio_num:.2f} / {ratio_den:.2f}")
    print(f"  L1/L2 = {lift_ratio:.4f}")
    print("\nThe final lift ratio L1/L2 as a fraction is 1/9.")

# Run the solver
if __name__ == "__main__":
    solve_aerofoil_problem()
    # To conform to the output format, we extract the numerical answer.
    # The calculated ratio is 1/9
    sys.stdout.write("<<<0.1111>>>")
