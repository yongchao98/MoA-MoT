def display_it3_vertical_slice_formula():
    """
    This function formulates and prints the mathematical equation for the
    upper bound of a vertical cross-section of an Interval Type-3 Membership Function (IT3 MF),
    based on a Gaussian paradigm with uncertain parameters.
    """

    # --- 1. Define Symbolic Representations ---
    # These strings represent the mathematical symbols in the equation.
    mu_upper_symbol = "μ_upper(u | x')"
    primary_variable_fixed = "x'"
    secondary_variable = "u"
    center_lower_func = f"c_lower({primary_variable_fixed})"
    center_upper_func = f"c_upper({primary_variable_fixed})"
    sigma_upper_func = f"σ_upper({primary_variable_fixed})"

    # --- 2. Construct the Piecewise Function ---
    # This represents the Upper Membership Function (UMF) of the vertical slice.
    # It is a Gaussian function with an uncertain mean, which creates a flat top.
    
    # Part 1: Left side of the Gaussian curve
    case1 = f"exp(-({secondary_variable} - {center_lower_func})^2 / (2 * {sigma_upper_func}^2))"
    condition1 = f"if {secondary_variable} < {center_lower_func}"

    # Part 2: Flat top of the function where membership is 1
    case2 = "1"
    condition2 = f"if {center_lower_func} <= {secondary_variable} <= {center_upper_func}"

    # Part 3: Right side of the Gaussian curve
    case3 = f"exp(-({secondary_variable} - {center_upper_func})^2 / (2 * {sigma_upper_func}^2))"
    condition3 = f"if {secondary_variable} > {center_upper_func}"

    # --- 3. Print Explanation of Each Term in the Equation ---
    print("The mathematical formulation for the upper bound of a vertical cross-section of an Interval Type-3 Membership Function (IT3 MF) is a piecewise function.")
    print("This function describes the Upper Membership Function (UMF) of the resulting Interval Type-2 fuzzy set for a fixed primary input x'.")
    print("\n--- Definition of Terms ---")
    print(f"  {mu_upper_symbol:<22}: The upper membership grade of the vertical slice.")
    print(f"  {secondary_variable:<22}: The secondary input variable, over which the cross-section is defined.")
    print(f"  {primary_variable_fixed:<22}: A fixed value of the primary input variable.")
    print(f"  {center_lower_func:<22}: The lower bound of the interval for the Gaussian's center at {primary_variable_fixed}.")
    print(f"  {center_upper_func:<22}: The upper bound of the interval for the Gaussian's center at {primary_variable_fixed}.")
    print(f"  {sigma_upper_func:<22}: The upper bound of the interval for the Gaussian's standard deviation at {primary_variable_fixed}, defining the maximum spread.")
    
    # --- 4. Print the Final Formatted Equation ---
    print("\n--- Final Equation ---")
    # We construct the final string with careful formatting for alignment.
    # The width for the first part is calculated to align the 'if' conditions.
    align_width = len(case1) + 2
    
    print(f"{mu_upper_symbol} = {{")
    print(f"  {case1 + ',':<{align_width}} {condition1}")
    print(f"  {case2 + ',':<{align_width}} {condition2}")
    print(f"  {case3 + ',':<{align_width}} {condition3}")
    print("}")

    # --- Prepare the single-line string for the final answer format ---
    final_answer_string = (
        f"μ_upper(u | x') = {{ "
        f"exp(-(u - c_lower(x'))^2 / (2 * σ_upper(x')^2)) if u < c_lower(x'); "
        f"1 if c_lower(x') <= u <= c_upper(x'); "
        f"exp(-(u - c_upper(x'))^2 / (2 * σ_upper(x')^2)) if u > c_upper(x') "
        f"}}"
    )
    return final_answer_string

# Execute the function to print the explanation and equation.
# The return value is captured for the final answer block.
final_answer = display_it3_vertical_slice_formula()
print(f"\n<<< {final_answer} >>>")