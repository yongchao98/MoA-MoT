import numpy as np

def calculate_hopf_charge():
    """
    Calculates the Hopf charge for the given vector field by determining
    the azimuthal and poloidal winding numbers.
    """
    # The Hopf charge Q_H for a field with this structure is the product of
    # two winding numbers: Q_H = n_f * n_G.

    # 1. Determine the azimuthal winding number n_f.
    # The field's azimuthal angle is defined as f = atan2(y, x).
    # By definition, atan2(y, x) increases by 2*PI for one full counter-clockwise
    # rotation around the z-axis. This corresponds to a winding number of 1.
    n_f = 1
    
    # 2. Determine the poloidal winding number n_G.
    # This number counts how many times the polar angle G covers the range [0, PI]
    # along a path from the axis of symmetry (z-axis) out to infinity.
    # We analyze this along a radial path in the z=0 plane, where rho = sqrt(x*x + y*y).
    #
    # The function G is defined as:
    # G = PI * exp(-10 * r2)
    # r2 = sqrt((rho*rho - 0.5)^2 + z*z)

    # At the center (rho=0, z=0):
    rho_center = 0.0
    r2_center = np.sqrt((rho_center**2 - 0.5)**2 + 0.0**2)
    G_center = np.pi * np.exp(-10 * r2_center)

    # On the hopfion core ring (where r2 is minimal, i.e., r2=0):
    # This occurs at rho=sqrt(0.5) and z=0.
    r2_core = 0.0
    G_core = np.pi * np.exp(-10 * r2_core)

    # At infinity (rho -> inf or z -> inf):
    # As r2 -> inf, exp(-10*r2) -> 0, so G -> 0.
    G_infinity = 0.0

    # As we travel from the center (rho=0) to infinity, the angle G starts
    # near 0 (G_center is very small), increases to a maximum of PI at the core,
    # and then decreases back to 0 at infinity. This trajectory sweeps the
    # range [0, PI] exactly once.
    # Therefore, the poloidal winding number n_G is 1.
    n_G = 1

    # 3. Calculate the total Hopf charge.
    Q_H = n_f * n_G

    print("The Hopf charge Q_H is the product of the azimuthal and poloidal winding numbers.")
    print(f"Azimuthal winding number n_f = {n_f}")
    print(f"Poloidal winding number n_G = {n_G}")
    print("\nThe final equation for the Hopf charge is:")
    print(f"{n_f} * {n_G} = {Q_H}")

calculate_hopf_charge()
<<<1>>>