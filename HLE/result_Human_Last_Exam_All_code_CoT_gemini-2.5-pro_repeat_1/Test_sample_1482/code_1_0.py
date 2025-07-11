import sympy as sp

def generate_it3_formulation():
    """
    Generates and prints the mathematical formulation for the UMF of a vertical
    cross-section of a Gaussian Interval Type-3 Membership Function.
    """
    # --- 1. Define Symbolic Variables ---
    # These symbols represent the variables and parameters in the formulation.
    u = sp.Symbol('u')  # Secondary input variable
    x_prime = sp.Symbol("x'")  # Fixed primary input variable
    c_u = sp.Symbol('c_u')  # Mean of the Gaussian function for the secondary variable u
    c_x = sp.Symbol('c_x')  # Mean of the Gaussian function for the primary variable x
    sigma_x_upper = sp.Symbol('sigma_x_upper', positive=True) # Std. dev. for the primary UMF
    sigma_u_min = sp.Symbol('sigma_u_min', positive=True) # Min std. dev. for the secondary variable u
    sigma_u_max = sp.Symbol('sigma_u_max', positive=True) # Max std. dev. for the secondary variable u

    # --- 2. Explain and Define Formulation Components ---
    print("This script provides the mathematical formulation for the Upper Membership Function (UMF)")
    print("of a vertical cross-section of a Gaussian Interval Type-3 Membership Function (IT3 MF).")
    print("-" * 70)
    print("A vertical cross-section at a fixed primary input x = x' is an Interval Type-2 (IT2) MF.")
    print("Its UMF, mu_upper(u | x'), is formulated as follows:\n")

    # Component 1: The UMF of the primary IT2 fuzzy set for variable x, evaluated at x'.
    # This function determines how much "blur" to apply to the secondary variable u.
    print("1. First, define the UMF of the primary IT2 set for x, evaluated at x'.")
    print("   This is a Gaussian function describing the upper bound of uncertainty for x.")
    mu_Ax_upper_expr = sp.exp(-sp.S(1)/2 * ((x_prime - c_x) / sigma_x_upper)**2)
    mu_Ax_upper_func = sp.Function('mu_A_x_upper')(x_prime)
    print("\nThe formula is:")
    sp.pprint(sp.Eq(mu_Ax_upper_func, mu_Ax_upper_expr), use_unicode=True)
    print("\n" + "-"*70 + "\n")

    # Component 2: The effective standard deviation for the vertical slice's UMF.
    # This value is linearly interpolated based on the result of Component 1.
    print("2. Next, calculate the effective standard deviation for the vertical slice's UMF.")
    print("   This standard deviation, sigma_upper(x'), is modulated by the primary membership value.")
    sigma_upper_expr = sigma_u_min + (sigma_u_max - sigma_u_min) * mu_Ax_upper_expr
    sigma_upper_func = sp.Function('sigma_upper')(x_prime)
    print("\nThe formula is:")
    sp.pprint(sp.Eq(sigma_upper_func, sigma_u_min + (sigma_u_max - sigma_u_min) * mu_Ax_upper_func), use_unicode=True)
    print("\n" + "-"*70 + "\n")

    # Component 3: The final UMF for the vertical cross-section.
    # This is a Gaussian function of 'u' using the standard deviation from Component 2.
    print("3. Finally, the UMF of the vertical cross-section is a Gaussian function of u.")
    print("   It uses the mean c_u and the effective standard deviation sigma_upper(x').")
    mu_upper_final_expr = sp.exp(-sp.S(1)/2 * ((u - c_u) / sigma_upper_expr)**2)
    mu_upper_func = sp.Function('mu_upper')(u, sp.S('x\'')) # Using a string for 'x|' notation
    print("\nThe final formulation is:")
    sp.pprint(sp.Eq(mu_upper_func, mu_upper_final_expr), use_unicode=True)
    print("\n" + "-"*70)
    print("In this equation, each part is defined as follows:")
    print(f"mu_upper(u, x') = {mu_upper_final_expr}")
    print(f"sigma_upper(x') = {sigma_upper_expr}")
    print(f"mu_A_x_upper(x') = {mu_Ax_upper_expr}")


if __name__ == '__main__':
    generate_it3_formulation()
