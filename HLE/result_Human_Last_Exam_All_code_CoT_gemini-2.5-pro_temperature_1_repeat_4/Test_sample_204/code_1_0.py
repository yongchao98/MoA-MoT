import numpy as np
from scipy.integrate import dblquad

def calculate_hopf_charge():
    """
    Calculates the Hopf charge of the given vector field by computing
    the linking number of the preimages of two points.
    
    The vector field n(x,y,z) is a map from R^3 to S^2.
    The Hopf charge H is the linking number Lk(n^{-1}(p1), n^{-1}(p2))
    for two regular values p1, p2 on S^2.

    1. Preimage of the South Pole n=(0,0,-1):
       cos(G) = -1 => G = PI => exp(-10*r2) = 1 => r2 = 0.
       This gives z=0 and x^2+y^2 = 0.5.
       This is a circle C_S with radius R_c = sqrt(0.5) in the xy-plane.

    2. Preimage of the equatorial point n=(1,0,0):
       G = PI/2 and f = 0.
       G = PI/2 => exp(-10*r2) = 0.5 => r2 = ln(2)/10.
       f = 0 => y=0, x>0.
       This gives (x^2-0.5)^2 + z^2 = (ln(2)/10)^2.
       This is a closed loop C_E in the xz-plane.

    3. The linking number is computed using the Gauss linking integral:
       Lk = (1 / (4*pi)) * integral_CS integral_CE [ (r1 - r2) . (dr1 x dr2) ] / |r1 - r2|^3 ds dt
    """

    # Constants from the problem definition
    R_c = np.sqrt(0.5)
    R_E = np.log(2) / 10.0

    # Parametrization of the curves C_S (r1) and C_E (r2)
    # C_S: r1(t) = (R_c*cos(t), R_c*sin(t), 0)
    # C_E: r2(s) = (sqrt(0.5 + R_E*cos(s)), 0, R_E*sin(s))

    def integrand(t, s):
        """
        Calculates the value of the integrand for the Gauss linking formula.
        
        Args:
            t: Parameter for curve C_S.
            s: Parameter for curve C_E.
        """
        # Position vectors
        r1 = np.array([R_c * np.cos(t), R_c * np.sin(t), 0])
        
        x2 = np.sqrt(0.5 + R_E * np.cos(s))
        z2 = R_E * np.sin(s)
        r2 = np.array([x2, 0, z2])
        
        # Derivatives of position vectors
        dr1 = np.array([-R_c * np.sin(t), R_c * np.cos(t), 0])
        
        # Avoid division by zero if x2 is somehow zero (it isn't for these params)
        if x2 < 1e-9:
            return 0
        dx2_ds = -R_E * np.sin(s) / (2 * x2)
        dz2_ds = R_E * np.cos(s)
        dr2 = np.array([dx2_ds, 0, dz2_ds])
        
        # Vector difference
        r_diff = r1 - r2
        
        # Denominator |r1 - r2|^3
        dist_cubed = np.sum(r_diff**2)**1.5
        if dist_cubed < 1e-9: # Should not happen as curves don't intersect
            return 0
            
        # Numerator (r1 - r2) . (dr1 x dr2) = det(r1-r2, dr1, dr2)
        numerator = np.dot(r_diff, np.cross(dr1, dr2))
        
        return numerator / dist_cubed

    # Perform the numerical double integration
    # The limits are 0 to 2*pi for both t and s
    integral_value, error = dblquad(integrand, 0, 2*np.pi, lambda t: 0, lambda t: 2*np.pi)
    
    # Calculate the linking number (Hopf Charge)
    factor = 1 / (4 * np.pi)
    hopf_charge = factor * integral_value

    print("The calculation follows the Whitehead formula for the Hopf charge,")
    print("which is computed here as the linking number of two preimages.")
    print("The linking number formula is: Lk = (1 / (4 * pi)) * Integral")
    print("\n--- Calculation Results ---")
    print(f"Integral Value: {integral_value}")
    print(f"Equation: {hopf_charge} = ({factor}) * ({integral_value})")
    print(f"The calculated Hopf charge is: {hopf_charge}")
    print(f"Rounding to the nearest integer, the Hopf charge is: {round(hopf_charge)}")

if __name__ == '__main__':
    calculate_hopf_charge()