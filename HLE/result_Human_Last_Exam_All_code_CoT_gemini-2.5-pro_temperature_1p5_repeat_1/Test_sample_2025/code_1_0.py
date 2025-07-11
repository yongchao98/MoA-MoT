import scipy.integrate as integrate
import numpy as np

def calculate_fz(z):
    """
    Calculates the value of the pdf f_Z(z) at a given point z.
    
    The calculation is based on the decomposition of the problem by the rank of X_1.
    f_Z(z) = 3*z*(1-z) + (1/8)*(g_12(z) + g_34(z) + 2*g_23(z)),
    where g_ij are pdfs of Z chosen from U[X_(i), X_(j)] from 4 points.
    """

    # Term 1: Contribution from X_1 being an outer point (X_(1) or X_(4)).
    term1 = 3 * z * (1 - z)

    # --- Numerical integration for the g_ij functions where j-i=1 ---
    
    # g_12(z) = integral from z to 1 of 12*(1-v)^2 * log(v/(v-z)) dv
    g12_integrand = lambda v: 12 * (1 - v)**2 * np.log(v / (v - z)) if v != z else 0
    g12_z, _ = integrate.quad(g12_integrand, z, 1)

    # g_34(z) = g_12(1-z) by symmetry.
    z_sym = 1 - z
    g34_integrand = lambda v: 12 * (1 - v)**2 * np.log(v / (v - z_sym)) if v != z_sym else 0
    g34_z, _ = integrate.quad(g34_integrand, z_sym, 1)
    
    # g_23(z) = integral from z to 1 of 24*(1-v)*(v*log(v/(v-z))-z) dv
    g23_integrand = lambda v: 24 * (1 - v) * (v * np.log(v / (v - z)) - z) if v != z else 0
    g23_z, _ = integrate.quad(g23_integrand, z, 1)

    # Term 2: Contribution from g_12(z)
    term2 = g12_z / 8
    
    # Term 3: Contribution from g_34(z)
    term3 = g34_z / 8
    
    # Term 4: Contribution from 2*g_23(z)
    term4 = (2 * g23_z) / 8

    # The final value is the sum of these terms.
    final_fz = term1 + term2 + term3 + term4

    # The problem asks to output each number in the final equation.
    print(f"The equation is: f(z) = 3*z*(1-z) + g_12(z)/8 + g_34(z)/8 + 2*g_23(z)/8")
    print(f"For z = {z}:")
    print(f"Term 1 (from outer points): 3*{z}*(1-{z}) = {term1:.8f}")
    print(f"Term 2 (from g_12({z})): {g12_z:.8f}/8 = {term2:.8f}")
    print(f"Term 3 (from g_34({z})): {g34_z:.8f}/8 = {term3:.8f}")
    print(f"Term 4 (from 2*g_23({z})): 2*{g23_z:.8f}/8 = {term4:.8f}")
    print(f"\nFinal Value f({z}) = {term1:.8f} + {term2:.8f} + {term3:.8f} + {term4:.8f} = {final_fz:.8f}")

if __name__ == '__main__':
    z_value = 0.2
    calculate_fz(z_value)
