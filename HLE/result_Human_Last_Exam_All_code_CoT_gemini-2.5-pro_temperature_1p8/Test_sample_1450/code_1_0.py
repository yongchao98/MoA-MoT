import numpy as np
from scipy import integrate
from scipy import special
from scipy import optimize

def fermi_dirac_integral(n, z):
    """
    Computes the complete Fermi-Dirac integral f_n(z).
    The definition is f_n(z) = 1/Gamma(n) * integral_0^inf t^(n-1)/(z^-1*exp(t) + 1) dt.
    """
    if not z > 0:
        return 0.0

    # The integrand function, written to be numerically stable.
    def integrand(t):
        log_z = np.log(z)
        val = t - log_z
        # For large (t - log_z), exp(val) can overflow.
        # In this case, the denominator is dominated by exp(val),
        # so the expression is approx. t^(n-1)*exp(-val), which is handled correctly by quad.
        # So no special clipping is needed here with modern scipy.
        return np.power(t, n - 1) / (np.exp(val) + 1.0)

    # Numerically integrate from 0 to infinity
    integral_val, _ = integrate.quad(integrand, 0, np.inf)

    # Normalize by the Gamma function Gamma(n)
    return integral_val / special.gamma(n)

def equation_to_solve(z):
    """
    This is the function for which we want to find the root.
    It represents f_{3/2}(z) - 0.75 * f_{5/2}(z) = 0.
    """
    if not z > 0:
        return np.nan # Return NaN for invalid fugacity values

    f_3_2 = fermi_dirac_integral(1.5, z)
    f_5_2 = fermi_dirac_integral(2.5, z)
    
    return f_3_2 - 0.75 * f_5_2

# --- Main execution ---
# Explain the problem and the derived equation
print("The problem requires finding the fugacity 'z' for a spin-1/2 ideal Fermi gas.")
print("The condition is that its number density is 75% of a classical ideal gas under the same pressure P and temperature T.")
print("This condition simplifies to the equation: f_{3/2}(z) = 0.75 * f_{5/2}(z)")
print("where f_n(z) is the Fermi-Dirac integral of order n.\n")
print("We solve the equation f_{3/2}(z) - 0.75 * f_{5/2}(z) = 0 numerically for z.\n")

# A large fugacity (z > 1) is expected, indicating quantum degeneracy.
# We bracket the root between 1 and 40 for the solver.
try:
    sol = optimize.root_scalar(equation_to_solve, bracket=[1, 40], method='brentq')
    if sol.converged:
        z_solution = sol.root

        # Verify the solution in the equation
        f_3_2_sol = fermi_dirac_integral(1.5, z_solution)
        f_5_2_sol = fermi_dirac_integral(2.5, z_solution)

        print(f"Solving for z...")
        print(f"Found solution at z = {z_solution:.4f}\n")
        
        print("Verifying the solution satisfies the equation f_{3/2}(z) = 0.75 * f_{5/2}(z):")
        print(f"Left side:  f_1.5({z_solution:.2f}) = {f_3_2_sol:.4f}")
        print(f"Right side: 0.75 * f_2.5({z_solution:.2f}) = {0.75 * f_5_2_sol:.4f}\n")
        
        # The problem asks for the answer to two significant digits.
        final_answer = f"{z_solution:.2g}"
        print(f"The value of the fugacity 'z' rounded to two significant digits is: {final_answer}")

    else:
        print("The root-finding algorithm did not converge.")

except ValueError as e:
    print(f"Could not solve the equation. The bracket might be incorrect. Error: {e}")

<<<28>>>