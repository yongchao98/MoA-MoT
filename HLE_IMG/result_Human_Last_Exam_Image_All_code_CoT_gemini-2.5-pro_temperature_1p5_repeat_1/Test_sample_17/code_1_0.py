import math

def demonstrate_stress_concentration():
    """
    This function demonstrates the concept of stress concentration at a sharp corner.
    """
    print("This script demonstrates the theoretical stress at a sharp corner (stress concentrator).")
    print("The maximum stress (sigma_max) at a notch is given by sigma_max = K_t * sigma_nominal.")
    print("Where K_t is the stress concentration factor and sigma_nominal is the applied stress.\n")
    print("For a sharp notch, K_t can be approximated by: K_t = 1 + 2 * sqrt(d / rho)")
    print("where 'd' is the notch depth and 'rho' is the radius of curvature at the notch tip.\n")

    # Let's assume some nominal values for demonstration purposes.
    # The final conclusion is independent of these specific values.
    sigma_nominal = 100  # Example nominal stress in MPa
    notch_depth_d = 2  # Example notch depth in mm

    print(f"Assuming a nominal applied stress of {sigma_nominal} MPa and a notch depth of {notch_depth_d} mm.")
    print("Let's see what happens to the maximum stress as the tip gets sharper (rho -> 0):\n")

    # A series of progressively smaller tip radii
    tip_radii_rho = [1.0, 0.1, 0.01, 0.001, 0.0001, 0.000001]

    print(f"{'Tip Radius (rho)':<20} | {'Stress Factor (K_t)':<25} | {'Max Stress (MPa)':<20}")
    print("-" * 70)

    for rho in tip_radii_rho:
        # Calculate the stress concentration factor
        K_t = 1 + 2 * math.sqrt(notch_depth_d / rho)
        
        # Calculate the maximum stress
        sigma_max = K_t * sigma_nominal
        
        print(f"{rho:<20.6f} | {K_t:<25.2f} | {sigma_max:<20.2f}")

    print("\n" + "-" * 70)
    print("As you can see, as the tip radius 'rho' approaches zero, the stress concentration factor 'K_t'")
    print("and the maximum stress 'sigma_max' increase without bound.\n")
    print("Therefore, for a theoretically perfectly sharp corner (rho = 0),")
    print("the theoretical stress at the tip is infinite.")


if __name__ == "__main__":
    demonstrate_stress_concentration()
