import math

def generate_it3_mf_formulation():
    """
    This script generates and explains the mathematical formulation for the upper bound
    of a vertical cross-section of a Gaussian Interval Type-3 Membership Function (IT3 MF).
    """

    # --- 1. Define the symbolic components of the equation ---
    lhs = "μ_upper(x, u)"
    func = "exp"
    secondary_var = "u"
    upper_bound_t2 = "μ_upper(x)"
    lower_bound_t2 = "μ_lower(x)"
    scale_param = "α_upper"

    # --- 2. Provide a clear explanation of the formulation and its terms ---
    print("This script provides the mathematical formulation for the upper bound of a vertical cross-section of a Gaussian Interval Type-3 Membership Function (IT3 MF).")
    print("This formulation describes the secondary membership grade at a fixed primary input 'x' as a function of the secondary variable 'u'.\n")

    print("--- Explanation of Terms ---")
    print(f"{lhs}: The final IT3 membership grade. It represents the upper bound of the vertical cross-section at a fixed 'x'.")
    print(f"{func}: The exponential function, forming the basis of the Gaussian shape.")
    print(f"{secondary_var}: The secondary input variable, representing a potential primary membership grade within the Footprint of Uncertainty (FOU).")
    print(f"{upper_bound_t2}: The upper membership function of the embedded Interval Type-2 set. For a Gaussian IT2 MF, this is exp(-((x-c)/σ_upper)^2).")
    print(f"{lower_bound_t2}: The lower membership function of the embedded Interval Type-2 set. For a Gaussian IT2 MF, this is exp(-((x-c)/σ_lower)^2).")
    print(f"{scale_param}: A user-defined scaling parameter (where {scale_param} ∈ (0, 1]) that controls the 'spread' or 'blurriness' of the secondary membership function.\n")

    # --- 3. Construct and print the full equation and its components ---
    numerator_inner = f"({secondary_var} - {upper_bound_t2})"
    denominator_inner = f"({scale_param} * ({upper_bound_t2} - {lower_bound_t2}))"
    exponent = f"-(({numerator_inner}) / ({denominator_inner}))^2"
    equation = f"{lhs} = {func}({exponent})"

    print("--- Mathematical Formulation ---")
    print(equation)
    print("\n--- Equation Components ---")
    print(f"Left-Hand Side (Result): {lhs}")
    print(f"Primary Function: {func}")
    print(f"Secondary Variable: {secondary_var}")
    print(f"Parameters from embedded Type-2 Set: {upper_bound_t2}, {lower_bound_t2}")
    print(f"Scaling Parameter: {scale_param}")
    print(f"Exponent Term: {exponent}")

    # --- 4. Output the final answer in the required format ---
    final_answer = f"<<<{equation}>>>"
    print(f"\n{final_answer}")

# Execute the function
generate_it3_mf_formulation()