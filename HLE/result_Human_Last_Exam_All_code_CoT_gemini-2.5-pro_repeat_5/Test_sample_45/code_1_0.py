import sympy as sp

def solve_telescope_problem():
    """
    This script derives the relationship between the focal length (f) of a
    liquid-mirror telescope and time (t) when driven by a constant power source.
    It then calculates the exponent 'n' in the relation f ∝ t^n.
    """

    # --- Step 1: Relate focal length (f) to angular velocity (ω) ---
    print("Step 1: Relate focal length (f) to angular velocity (ω).")
    print("The surface of a rotating fluid under gravity forms a paraboloid. The equation for the surface height z as a function of radius r is:")
    print("z(r) = (ω^2 * r^2) / (2 * g) + z_0")
    print("where g is the acceleration due to gravity and z_0 is the height at the center.\n")

    print("This is a parabola of the form y = a*x^2. The focal length 'f' of such a parabola is given by f = 1 / (4*a).")
    print("By comparing the forms, we identify the coefficient 'a' as: a = ω^2 / (2*g).")
    print("Substituting 'a' into the focal length formula gives:")
    print("f = 1 / (4 * (ω^2 / (2*g))) = g / (2 * ω^2)")
    print("Since g is a constant, the focal length is inversely proportional to the square of the angular velocity:")
    print("f ∝ 1 / ω^2  or  f ∝ (ω^2)^(-1)\n")

    # --- Step 2: Relate angular velocity (ω) to time (t) ---
    print("Step 2: Relate angular velocity (ω) to time (t).")
    print("The system is driven by a constant power source, P.")
    print("Power (P) is the rate of change of rotational kinetic energy (K): P = dK/dt.")
    print("The kinetic energy is K = (1/2) * I * ω^2, where I is the constant moment of inertia.\n")

    print("Substituting K into the power equation:")
    print("P = d/dt [ (1/2) * I * ω^2 ] = (1/2) * I * d(ω^2)/dt")
    print("This is a differential equation. We solve it by integrating with respect to time, starting from rest (ω=0 at t=0):")
    print("∫(from 0 to t) P dt' = ∫(from 0 to K) dK'")
    print("P * t = K = (1/2) * I * ω^2\n")

    print("Rearranging for ω^2, we get:")
    print("ω^2 = (2 * P * t) / I")
    print("Since P and I are constants, the square of the angular velocity is directly proportional to time:")
    print("ω^2 ∝ t\n")

    # --- Step 3: Combine relationships to find the exponent n ---
    print("Step 3: Combine the relationships to find the final dependency of f on t.")
    print("From Step 1, we have the proportionality: f ∝ (ω^2)^(-1)")
    print("From Step 2, we have the proportionality: ω^2 ∝ t\n")

    print("Substituting the time dependency of ω^2 into the focal length equation:")
    print("f ∝ (t)^(-1)")
    print("This can be written as: f ∝ t^n, where n is the exponent we need to find.\n")

    # Final Calculation
    n = -1
    print("By comparing the final expression f ∝ t^(-1) with f ∝ t^n, we find the value of n.")
    print(f"The final equation is f ∝ t^({n})")
    print(f"Therefore, the value of the exponent n is: {n}")

if __name__ == '__main__':
    solve_telescope_problem()