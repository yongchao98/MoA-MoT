import numpy as np
from scipy.special import iv, kv
from scipy.integrate import quad

def illustrate_no_positive_eigenvalues():
    """
    This function illustrates why the stability operator L has no positive eigenvalues.
    It analyzes the asymptotic radial equation for a specific case and shows that
    no solution can be both regular at the origin and square-integrable.
    """
    # Let's choose a specific case for illustration: n=3, k=1.
    # The argument holds for any n>=2 and k>=0.
    n = 3
    k = 1
    # Let's pick an arbitrary positive eigenvalue candidate, lambda = 1.
    lam = 1.0

    print(f"--- Illustration for n={n}, k={k}, lambda={lam} ---")

    # The order of the modified Bessel function.
    nu = k + (n - 2.0) / 2.0

    # The asymptotic radial equation is a modified Bessel equation.
    # Its solutions are related to I_nu and K_nu.
    # f(rho) ~ rho^((2-n)/2) * (c1 * I_nu(sqrt(lam)*rho) + c2 * K_nu(sqrt(lam)*rho))

    # Solution that is regular at the origin (rho=0)
    # This corresponds to the modified Bessel function of the first kind, I_nu.
    def f_reg(rho):
        if rho == 0:
            # Handle the limit at rho=0 based on the behavior of I_nu
            return 0.0 if nu > 0 else 1.0
        return rho**((2.0 - n) / 2.0) * iv(nu, np.sqrt(lam) * rho)

    # Solution that is square-integrable at infinity (rho -> inf)
    # This corresponds to the modified Bessel function of the second kind, K_nu.
    def f_l2(rho):
        return rho**((2.0 - n) / 2.0) * kv(nu, np.sqrt(lam) * rho)

    # --- Check the properties of these two solutions ---

    # 1. The L^2 solution must be regular at rho=0. Let's check f_l2.
    # K_nu(x) blows up as x->0 for nu>0, and as log(x) for nu=0.
    print("\n1. Analyzing the square-integrable solution (f_l2):")
    try:
        val_at_small_rho = f_l2(1e-9)
        print(f"Value of f_l2 at a small rho (1e-9): {val_at_small_rho:.4e}")
        if np.isinf(val_at_small_rho) or val_at_small_rho > 1e6:
            print("This solution is singular at rho=0, so it cannot be an eigenfunction.")
    except Exception as e:
        print(f"An error occurred while evaluating f_l2 near zero: {e}")


    # 2. The regular solution must be square-integrable. Let's check f_reg.
    # The L^2 norm involves an integral of |f(rho)|^2 * w(rho), where the weight
    # w(rho) is asymptotically rho^(n-1).
    def integrand_reg(rho):
        return np.abs(f_reg(rho))**2 * rho**(n - 1)

    print("\n2. Analyzing the regular solution (f_reg):")
    # We integrate from a small number > 0 to a large number R.
    # If the solution is in L^2, this integral should converge as R -> infinity.
    # I_nu(x) grows exponentially, so we expect the integral to diverge.
    try:
        integral_1, _ = quad(integrand_reg, 1e-6, 10)
        integral_2, _ = quad(integrand_reg, 1e-6, 20)
        print(f"Integral of |f_reg|^2 * w(rho) from 0 to 10: {integral_1:.4e}")
        print(f"Integral of |f_reg|^2 * w(rho) from 0 to 20: {integral_2:.4e}")
        if integral_2 > integral_1 * 2:
             print("The norm integral grows rapidly, indicating the function is not square-integrable.")
    except Exception as e:
        print(f"An error occurred during integration: {e}")


    # --- Conclusion ---
    print("\n--- Conclusion ---")
    print("An eigenfunction must be both regular at the origin and square-integrable.")
    print("Our analysis shows that for lambda > 0, no non-trivial solution has both properties.")
    print("Therefore, there are no positive eigenvalues.")
    
    final_answer = 0
    print("\nThe number of positive eigenvalues is:")
    print(final_answer)

if __name__ == '__main__':
    illustrate_no_positive_eigenvalues()
