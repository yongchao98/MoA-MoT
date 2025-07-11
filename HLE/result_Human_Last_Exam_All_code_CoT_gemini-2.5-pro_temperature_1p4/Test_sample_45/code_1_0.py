def solve_focal_length_exponent():
    """
    This function determines the exponent 'n' in the relationship f ∝ t^n
    for a liquid-mirror telescope rotated by a constant power source.

    The derivation steps are as follows:
    1. The focal length (f) of a rotating liquid mirror is related to its angular
       speed (ω) by the equation: f = g / (2 * ω^2).
       This implies: f ∝ 1 / ω^2

    2. A constant power source (P) provides energy to the system. Power is the
       rate of change of kinetic energy (K): P = dK/dt.
       Rotational kinetic energy is K = (1/2) * I * ω^2, where I is the
       moment of inertia.

    3. Since P and I are constant, we have d(ω^2)/dt = constant.
       Integrating from rest (ω=0 at t=0), we find that ω^2 is proportional to t.
       This implies: ω^2 ∝ t

    4. Combining the two proportionalities:
       f ∝ 1 / ω^2  and  ω^2 ∝ t
       Therefore, f ∝ 1 / t, which can be written as f ∝ t^(-1).

    5. By comparing f ∝ t^n with f ∝ t^(-1), we find the value of n.
    """
    
    # The final relation is f ∝ t^n, and we found f ∝ t^(-1).
    n = -1
    
    print("The relationship between focal length (f) and time (t) is f ∝ t^n.")
    print("Based on the physics of a constant power source, we derived that f ∝ t^(-1).")
    print("Therefore, the value of the exponent n is:")
    
    # Output the final number in the equation, which is the exponent n.
    print(n)

solve_focal_length_exponent()
<<<-1>>>