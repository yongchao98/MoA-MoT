import math

def solve_mass_ratio():
    """
    Calculates the leading-order asymptotic mass ratio for CP(N-1) models.
    """
    # The problem asks for the leading-order asymptotic mass ratio between the lightest solitonic excitation
    # and the subsequent higher excitation in the spectrum of CP(N-1) models as N approaches infinity.

    # Step 1: The mass spectrum of these excitations is given by m_k = m_1 * sin(k*pi/N) / sin(pi/N).
    # The lightest excitation is the k=1 state (mass m_1).
    # The subsequent higher excitation is the k=2 state (mass m_2).

    # Step 2: We need to find the ratio m_2 / m_1.
    # Ratio = [m_1 * sin(2*pi/N) / sin(pi/N)] / m_1 = sin(2*pi/N) / sin(pi/N)

    # Step 3: This expression simplifies using the identity sin(2x) = 2*sin(x)*cos(x)
    # to 2*cos(pi/N).

    print("The mass ratio R between the second and first excitation is given by the formula:")
    print("R = 2 * cos(pi/N)")
    print("\nWe need to find the limit of this ratio as N approaches infinity.")

    # Step 4: We calculate the limit: lim_{N->inf} [2 * cos(pi/N)].
    # As N -> inf, the term pi/N approaches 0.
    limit_of_pi_over_N = 0
    cos_of_limit = math.cos(limit_of_pi_over_N)
    
    # The factor '2' remains a constant.
    factor = 2

    final_ratio = factor * cos_of_limit

    print("\nThe final equation is derived from the limit calculation:")
    print(f"{factor} * cos({limit_of_pi_over_N}) = {factor} * {int(cos_of_limit)} = {int(final_ratio)}")

    print("\nThe leading-order asymptotic mass ratio is:")
    print(int(final_ratio))

solve_mass_ratio()
<<<2>>>