import math

def demonstrate_unboundedness():
    """
    This function demonstrates the concept from question (a).
    It shows that for a supercritical exponent p, the functional J_t
    becomes unbounded from below as t -> +infinity.
    """
    # Let's choose a value for the anisotropy parameter s, assuming s >= 1.
    s = 2.0

    # Calculate the critical exponent for p based on the derived formula.
    p_crit = 2 * (3 * s + 1) / (s + 1)

    # Choose a value for p that is in the supercritical range (p > p_crit).
    p = 5.0

    print(f"Chosen parameters: s = {s}, p = {p}")
    print(f"The critical exponent is p_crit = 2*(3*s + 1)/(s + 1) = {p_crit:.4f}")

    if p > p_crit:
        print(f"The condition p > p_crit ({p} > {p_crit:.4f}) is satisfied.")
    else:
        print(f"The condition p > p_crit ({p} > {p_crit:.4f}) is NOT satisfied.")
        return

    # For t -> infinity, with s > 1, the scaled functional J(u_t, v_t) behaves like:
    # J(t) approx C1 * t^(2s) - C2 * t^((s+1)*(p/2-1))
    # where C1 and C2 are positive constants related to the norms of u and v.
    # Let's model this with C1=1, C2=1 for simplicity.
    
    gamma_kinetic = 2 * s
    gamma_potential = (s + 1) * (p / 2 - 1)

    def J_scaled(t):
        return t**gamma_kinetic - t**gamma_potential

    print(f"\nModeling the scaled functional as J(t) = t^({gamma_kinetic}) - t^({gamma_potential})")
    
    # Check the value for a large t
    t_large = 10.0
    j_val = J_scaled(t_large)

    print(f"For a large t = {t_large}, the value is J({t_large}) = {j_val:.2f}")
    if j_val < 0:
        print("As t -> +infinity, the functional becomes unbounded from below.")

    # As requested, output each number in the final equation p > 2*(3*s+1)/(s+1)
    print("\nThe numbers in the inequality p > 2 * (3 * s + 1) / (s + 1) are:")
    # The equation is 5 > 2 * (3 * 2 + 1) / (2 + 1)
    numbers = [p, 2, 3, s, 1, s, 1]
    # To avoid floats in the final list for clarity, casting s to int for this print
    numbers_str = [str(int(n)) if n == float(int(n)) else str(n) for n in numbers]
    print(', '.join(numbers_str))


demonstrate_unboundedness()