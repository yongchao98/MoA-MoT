import math

def solve_raindrop_fall():
    """
    This function calculates the time for the raindrop to fall,
    following the physical derivation.

    The derivation steps are:
    1. The raindrop's mass (m) and rate of mass gain (dm/dt) are modeled:
       m = ρ * (4/3) * π * r^3
       dm/dt = Ρ * π * r^2 * v
    2. Newton's second law for variable mass is F_ext = d(mv)/dt.
       With F_ext = mg, this gives mg = m*a + v*(dm/dt).
    3. Assuming the raindrop starts with negligible size (r_0 = 0), these
       equations can be solved to find that the acceleration 'a' is constant.
       The densities ρ and Ρ cancel out, yielding:
       a = g / 7
    4. The time 't' to fall a height 'H' from rest with constant acceleration 'a'
       is given by the kinematic equation: H = (1/2)*a*t^2.
       t = sqrt(2 * H / a)
    5. Substituting a = g/7, we get:
       t = sqrt(2 * H / (g/7)) = sqrt(14 * H / g)
    6. The initial height is interpreted as H = 350 m.
    """

    # Initial height in meters
    H = 350

    # The numerator in the expression for t^2 = (14 * H) / g
    numerator_val = 14 * H

    # The square root of the numerator
    final_numerator = math.sqrt(numerator_val)

    # Print the final equation and its simplified form
    print("The final equation for the time of fall 't' is derived from t = sqrt(14 * H / g).")
    print(f"Substituting the height H = {H} m:")
    print(f"t = sqrt(14 * {H} / g)")
    print(f"t = sqrt({numerator_val} / g)")
    print(f"t = {int(final_numerator)} / sqrt(g)")

solve_raindrop_fall()