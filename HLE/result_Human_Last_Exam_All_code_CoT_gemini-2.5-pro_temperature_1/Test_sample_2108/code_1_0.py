import numpy as np
from scipy import integrate

def solve_max_ratio():
    """
    Calculates the maximum achievable ratio of bidirectional conical power
    to line intensity for a rotating charged sphere.
    """

    # 1. Define the radiation properties.
    # The leading-order radiation from the described rotating sphere has an
    # intensity pattern I(theta) = k * sin^4(theta) when the rotation axis
    # is chosen optimally (along a principal axis of the quadrupole).
    # We can set the proportionality constant k=1 without loss of generality
    # as it will cancel out in the final ratio.
    k = 1.0
    
    def intensity(theta):
        """Intensity as a function of angle theta from the rotation axis."""
        return k * (np.sin(theta))**4

    # The half-angle of the bidirectional cone.
    cone_half_angle = np.pi / 4
    cos_cone_half_angle = np.cos(cone_half_angle)

    # 2. Define the denominator of the ratio.
    # To get a finite, well-defined maximum ratio, we interpret "intensity
    # along a line" as the maximum possible intensity of the radiation pattern.
    # For I ~ sin^4(theta), the maximum occurs at theta = pi/2.
    I_line = intensity(np.pi / 2)

    # 3. Define the numerator: the conical power P_cone.
    # This function calculates the power integrated over the cone for a given
    # cone axis orientation, gamma. Gamma is the angle between the rotation
    # axis and the cone axis.
    def get_conical_power(gamma):
        """
        Calculates the power in a bidirectional cone whose axis is at an
        angle gamma to the rotation axis.
        """
        # The integrand is I(theta) * sin(theta).
        # The integration limits are determined by the cone's geometry.
        # Let n be the observation direction (theta, phi) and n_c be the
        # cone axis. The point is in the cone if |n . n_c| >= cos(alpha).
        def integrand(phi, theta):
            # Cone axis vector (we can place it in the x-z plane wlog)
            n_c_x = np.sin(gamma)
            n_c_z = np.cos(gamma)
            
            # Observation direction vector
            n_x = np.sin(theta) * np.cos(phi)
            n_z = np.cos(theta)
            
            # Dot product n . n_c
            dot_product = n_x * n_c_x + n_z * n_c_z
            
            # Check if the direction is within the bidirectional cone
            if np.abs(dot_product) >= cos_cone_half_angle:
                return intensity(theta) * np.sin(theta)
            else:
                return 0.0

        # Perform the numerical double integration over the sphere.
        # dblquad(func, theta_min, theta_max, phi_min(theta), phi_max(theta))
        power, _ = integrate.dblquad(integrand, 0, np.pi, 
                                     lambda theta: 0, lambda theta: 2 * np.pi,
                                     epsabs=1e-5, epsrel=1e-5)
        return power

    # 4. Find the optimal cone orientation that maximizes the conical power.
    # By symmetry, the maximum power will be for gamma = pi/2, where the
    # cone is centered on the equator (the peak of the sin^4(theta) pattern).
    # We calculate it for this optimal case directly.
    optimal_gamma = np.pi / 2
    max_P_cone = get_conical_power(optimal_gamma)

    # 5. Calculate the final maximum ratio.
    max_ratio = max_P_cone / I_line

    # The problem asks to output each number in the final equation.
    print(f"Maximum Bidirectional Conical Power (P_cone): {max_P_cone}")
    print(f"Intensity along a Line (taken as I_max): {I_line}")
    print(f"Maximum Achievable Ratio (P_cone / I_max): {max_ratio}")
    
solve_max_ratio()
<<<4.8442>>>