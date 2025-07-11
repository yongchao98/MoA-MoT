import numpy as np
from scipy.special import jn_zeros

def solve_critical_chemical_potential():
    """
    Calculates the critical chemical potential for scalar condensation in a hard-wall
    holographic model with an Einstein-Gauss-Bonnet background.

    The calculation is based on the following assumptions:
    - The model is a 5D hard-wall model at zero temperature.
    - The IR cutoff is set to z_IR = 1, which defines the energy scale.
    - The Gauss-Bonnet coupling lambda_GB = 0.1.
    - The scalar field has mass squared m^2 = -3 (dual to a Delta=3 operator).
    - The scalar field has charge q = 1.
    """
    # Parameters from the problem and standard assumptions
    lambda_gb = 0.1
    m_squared = -3.0
    q = 1.0
    z_ir = 1.0  # IR cutoff, setting the energy scale

    # At T=0, the metric function f(z) becomes a constant f_0
    f0 = (1 - np.sqrt(1 - 4 * lambda_gb)) / (2 * lambda_gb)

    # The scalar field EOM at T=0 simplifies to a Bessel equation.
    # The order of the Bessel function, nu, depends on m^2 and f0.
    # The EOM is psi'' + (f'/f - 3/z)psi' + (q^2*phi^2/f^2 - m^2/(z^2*f))psi = 0.
    # For T=0, f->f0, phi->mu. A change of variables psi=z^2*chi gives a Bessel equation for chi.
    # The order of the Bessel function is nu = sqrt(4 + m^2/f0).
    nu_squared = 4 + m_squared / f0
    if nu_squared < 0:
        print("Error: nu^2 is negative. Condensation might not occur in this model.")
        return

    nu = np.sqrt(nu_squared)

    # The hard-wall boundary condition psi(z_IR) = 0 implies that the argument
    # of the Bessel function must be a zero of J_nu. We take the first zero
    # for the critical (lowest) chemical potential.
    # The argument is k*z_IR, where k = q*mu/f0.
    # So, (q * mu_c / f0) * z_IR = j_{nu,1}, where j_{nu,1} is the first zero of J_nu(x).
    
    # Find the first zero of the Bessel function of order nu
    num_zeros = 1
    first_zero = jn_zeros(nu, num_zeros)[0]

    # Calculate the critical chemical potential mu_c
    mu_c = (f0 * first_zero) / (q * z_ir)
    
    # Print the step-by-step calculation
    print(f"1. Set parameters:")
    print(f"   Gauss-Bonnet coupling lambda_GB = {lambda_gb}")
    print(f"   Scalar mass squared m^2 = {m_squared}")
    print(f"   Scalar charge q = {q}")
    print(f"   IR cutoff z_IR = {z_ir} (sets the energy scale)")
    print("")

    print(f"2. Calculate the T=0 metric factor f0:")
    print(f"   f0 = (1 - sqrt(1 - 4 * {lambda_gb})) / (2 * {lambda_gb})")
    print(f"   f0 = {f0:.4f}")
    print("")

    print(f"3. Determine the order nu of the Bessel function:")
    print(f"   nu^2 = 4 + m^2 / f0 = 4 + {m_squared} / {f0:.4f} = {nu_squared:.4f}")
    print(f"   nu = sqrt({nu_squared:.4f}) = {nu:.4f}")
    print("")

    print(f"4. Find the first zero of the Bessel function J_nu(x):")
    print(f"   First zero j_{{{nu:.3f}},1}} = {first_zero:.4f}")
    print("")

    print(f"5. Calculate the critical chemical potential mu_c:")
    print(f"   mu_c = (f0 * first_zero) / (q * z_IR)")
    print(f"   mu_c = ({f0:.4f} * {first_zero:.4f}) / ({q} * {z_ir})")
    print(f"   mu_c = {mu_c:.4f}")
    print("")
    
    # Final equation with all numbers
    final_equation = f"mu_c = ((1 - sqrt(1 - 4 * {lambda_gb})) / (2 * {lambda_gb})) * {first_zero:.4f} / ({q} * {z_ir})"
    final_value_calc = f"mu_c = {f0:.4f} * {first_zero:.4f} = {mu_c:.4f}"
    
    print("Final equation with numbers plugged in:")
    print(f"{mu_c:.4f} = ((1 - sqrt(1 - 4 * {lambda_gb})) / (2 * {lambda_gb})) * {first_zero:.4f}")


solve_critical_chemical_potential()