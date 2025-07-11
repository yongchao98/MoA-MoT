import math

def calculate_hopf_charge():
    """
    Calculates the Hopf charge for the given vector field based on its topological structure.
    The Hopf Charge H is the product of the toroidal (k) and poloidal (P) winding numbers.
    H = k * P
    """

    # 1. Determine the toroidal winding number, k.
    # The vector field uses f = atan2(y,x), which wraps once around the z-axis.
    # This corresponds to a toroidal winding number k = 1.
    k = 1

    # 2. Determine the poloidal winding number, P.
    # P = (G_at_infinity - G_at_core) / PI
    # G = PI * exp(-10 * r2)

    # At spatial infinity, r2 -> infinity, so G -> 0.
    G_at_infinity = 0

    # At the core of the hopfion, r2 -> 0, so G -> PI * exp(0) = PI.
    G_at_core = math.pi
    
    # Calculate P. The result must be an integer.
    P = (G_at_infinity - G_at_core) / math.pi
    P_int = int(P)

    # 3. Calculate the total Hopf charge, H.
    H = k * P_int
    H_int = int(H)
    
    # Print the breakdown of the final calculation as requested.
    print("The Hopf charge is calculated as H = k * P")
    print(f"Toroidal winding number k = {k}")
    print(f"Poloidal winding number P = {P_int}")
    print(f"Final equation: H = {k} * {P_int} = {H_int}")


if __name__ == "__main__":
    calculate_hopf_charge()