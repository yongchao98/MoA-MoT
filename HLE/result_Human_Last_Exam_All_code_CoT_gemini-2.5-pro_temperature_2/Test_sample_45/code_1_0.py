import sympy

def solve_telescope_problem():
    """
    Derives the exponent n in the relationship f ∝ t^n for a liquid-mirror telescope.

    This function follows these physical principles:
    1.  Relates the focal length (f) to the angular speed (ω) of the rotating liquid.
    2.  Relates the angular speed (ω) to time (t) under constant power input.
    3.  Combines these relationships to find the final dependency of f on t.
    """

    print("Step 1: Find the relationship between focal length (f) and angular speed (ω)")
    # The surface of a rotating liquid forms a paraboloid with the equation:
    # z = (ω^2 * r^2) / (2 * g)
    # The equation for a parabolic mirror is r^2 = 4 * f * z.
    # By equating r^2/z from both equations, we get:
    # 4 * f = (2 * g) / ω^2
    # f = g / (2 * ω^2)
    # This means f is proportional to ω^(-2).
    n_f_vs_omega = -2
    print("The shape of the rotating liquid is a parabola, leading to the relationship:")
    print(f"f ∝ ω^a, where the exponent a = {n_f_vs_omega}\n")

    print("Step 2: Find the relationship between angular speed (ω) and time (t)")
    # A constant power source (P) provides energy to the system.
    # The work done is W = P * t. This becomes the rotational kinetic energy (K).
    # K = P * t
    # The kinetic energy of a rotating body is K = (1/2) * I * ω^2.
    # Equating the two expressions for K:
    # (1/2) * I * ω^2 = P * t
    # ω^2 ∝ t^1
    # This means ω is proportional to t^(1/2).
    n_omega_vs_t = 1/2
    print("The constant power source leads to a relationship for angular speed over time:")
    print("ω^2 ∝ t^b, where the exponent b = 1")
    print(f"This can be expressed as: ω ∝ t^c, where the exponent c = {n_omega_vs_t}\n")

    print("Step 3: Combine the relationships to find f ∝ t^n")
    # We have f ∝ ω^(-2) and ω ∝ t^(1/2).
    # Substitute the expression for ω into the expression for f:
    # f ∝ (t^(1/2))^(-2)
    # f ∝ t^(1/2 * -2)
    # f ∝ t^(-1)
    # So, the final exponent n is -1.
    final_n = n_f_vs_omega * n_omega_vs_t
    print("To find the final exponent n for the relation f ∝ t^n, we combine the exponents:")
    print("n = (exponent a) * (exponent c)")
    print(f"The equation with the numbers is: n = ({n_f_vs_omega}) * ({n_omega_vs_t})")
    print(f"The final calculated exponent n = {int(final_n)}\n")

    print(f"Therefore, the focal length f is proportional to t^{int(final_n)}.")

solve_telescope_problem()
print("<<<"+str(int(-1))+">>>")