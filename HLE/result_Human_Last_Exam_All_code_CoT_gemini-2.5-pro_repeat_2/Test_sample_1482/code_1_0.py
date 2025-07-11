import math

def generate_it3_formulation():
    """
    This function generates and prints the mathematical formulation for the vertical
    cross-section of a Gaussian-based Interval Type-3 Membership Function (IT3 MF).
    """
    # 1. Define the parameters for the Gaussian functions.
    # These parameters define the primary uncertainty of the IT3 MF.
    c = 5.0
    sigma_upper = 1.0  # Std dev for the upper bound of the primary FOU
    sigma_lower = 2.0  # Std dev for the lower bound of the primary FOU
    
    # This parameter controls the spread (uncertainty) of the vertical cross-section.
    k = 0.25
    
    # 2. Print the parameters for clarity.
    print("Parameters for the Gaussian IT3 MF formulation:")
    print(f"Center (c): {c}")
    print(f"Upper Bound Standard Deviation (σ_upper): {sigma_upper}")
    print(f"Lower Bound Standard Deviation (σ_lower): {sigma_lower}")
    print(f"Spread Control Parameter (k): {k}")
    print("-" * 60)

    # 3. Define the components of the final equation as formatted strings.
    # These represent the Gaussian functions for the primary FOU bounds.
    mu_upper_x_str = f"exp(-0.5 * ((x - {c}) / {sigma_upper})**2)"
    mu_lower_x_str = f"exp(-0.5 * ((x - {c}) / {sigma_lower})**2)"

    # These define the center and spread of the Gaussian for the vertical cross-section.
    center_Cx_str = f"({mu_upper_x_str} + {mu_lower_x_str}) / 2.0"
    spread_Sx_str = f"{k} * ({mu_upper_x_str} - {mu_lower_x_str})"

    # 4. Construct and print the full mathematical formulation.
    # This represents the upper bound of the FOU of the vertical cross-section.
    print("Mathematical Formulation for the Upper Bound of the Vertical Cross-Section μ_Ã_upper(x, u):")
    print("\nμ_Ã_upper(x, u) = exp( -0.5 * ( (u - C(x)) / S(x) )**2 )")

    print("\nWhere the components are defined as:")
    print(f"\nC(x) = {center_Cx_str}")
    print(f"\nS(x) = {spread_Sx_str}")
    print("-" * 60)

    # 5. Extract and print all numbers used in the final combined equation for clarity.
    print("Individual numerical constants in the full formulation:")
    print("From the Gaussian exponential function: -0.5")
    print("From the squaring operation: 2")
    print(f"From the center C(x) definition: 2.0")
    print(f"From the center c of the primary Gaussians: {c}")
    print(f"From the standard deviation σ_upper: {sigma_upper}")
    print(f"From the standard deviation σ_lower: {sigma_lower}")
    print(f"From the spread control parameter k: {k}")

# Execute the function to print the formulation.
if __name__ == "__main__":
    generate_it3_formulation()