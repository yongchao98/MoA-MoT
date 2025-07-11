def solve_liquid_mirror_telescope():
    """
    This function calculates the exponent 'n' in the relationship f ∝ t^n
    for a liquid-mirror telescope driven by a constant power source.

    The derivation is as follows:

    1. The focal length (f) of a rotating liquid mirror is related to its angular velocity (ω) by:
       f = g / (2 * ω^2), where g is the acceleration due to gravity.
       This implies: f ∝ 1 / ω^2  or  f ∝ ω^(-2)

    2. The constant power source (P) relates to the angular velocity (ω) and moment of inertia (I) by:
       P = τ * ω = (I * α) * ω = I * (dω/dt) * ω
       Separating variables and integrating from t=0 (where ω=0):
       P * t = (1/2) * I * ω^2
       This implies: ω^2 ∝ t

    3. Combining the two proportionalities:
       f ∝ (ω^2)^(-1)
       f ∝ t^(-1)

    4. Comparing this with f ∝ t^n, we find that n = -1.
    """
    print("The relationship between the focal length f and time t is derived as f ∝ t^n.")
    print("We need to find the value of the exponent n.")

    # The final equation is f ∝ t^(-1)
    n = -1

    print("\nThe final equation has the form f ∝ t^n, where the exponent is:")
    print(n)

solve_liquid_mirror_telescope()
