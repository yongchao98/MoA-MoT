import math

def generate_it3_mf_formulation():
    """
    Generates and prints the mathematical formulation for the vertical
    cross-sections of a Gaussian-based Interval Type-3 Membership Function.
    """

    # --- Define Numerical Parameters for the Example ---

    # The fixed value for the primary input variable x for the vertical cross-section.
    x_val = 6.0

    # The center (mean) of the Gaussian functions.
    c = 5.0

    # Parameters for the linear model of the lower bound standard deviation: sigma_L(u) = s_L1*u + s_L0
    s_L1 = 0.5
    s_L0 = 1.0

    # Parameters for the linear model of the upper bound standard deviation: sigma_U(u) = s_U1*u + s_U0
    # We choose parameters such that sigma_U(u) >= sigma_L(u) for u in [0, 1].
    s_U1 = 0.5
    s_U0 = 1.5

    # --- Construct and Print the Formulations ---

    print("The mathematical formulation for the uncertainty bounds of the vertical cross-section is as follows:\n")

    # Formulation for the Upper Bound (using the smaller standard deviation, sigma_L)
    # The term (x_val - c)^2 is pre-calculated for clarity in the final equation.
    numerator_val = (x_val - c)**2
    
    print("Upper Bound Membership Function:")
    # The final equation shows each number explicitly.
    upper_bound_eq = f"mu_bar(x={x_val}, u) = exp(-0.5 * (({x_val} - {c}) / ({s_L1}*u + {s_L0}))^2)"
    upper_bound_simplified = f"                  = exp(-{0.5 * numerator_val} / ({s_L1}*u + {s_L0})^2)"
    print(upper_bound_eq)
    print(upper_bound_simplified)

    print("\n" + "="*80 + "\n")

    # Formulation for the Lower Bound (using the larger standard deviation, sigma_U)
    print("Lower Bound Membership Function:")
    # The final equation shows each number explicitly.
    lower_bound_eq = f"mu_underline(x={x_val}, u) = exp(-0.5 * (({x_val} - {c}) / ({s_U1}*u + {s_U0}))^2)"
    lower_bound_simplified = f"                      = exp(-{0.5 * numerator_val} / ({s_U1}*u + {s_U0})^2)"
    print(lower_bound_eq)
    print(lower_bound_simplified)


if __name__ == "__main__":
    generate_it3_mf_formulation()
