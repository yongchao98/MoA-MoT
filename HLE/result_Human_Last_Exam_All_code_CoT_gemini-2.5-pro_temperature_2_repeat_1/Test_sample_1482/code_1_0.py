import sys

def formulate_it3_vertical_slice():
    """
    This function generates and explains the mathematical formulation for a vertical
    cross-section of an Interval Type-3 Membership Function using a Gaussian paradigm.
    """

    # --- 1. Introduction ---
    print("In an Interval Type-3 (IT3) Fuzzy Logic System, a vertical cross-section is examined by fixing the primary input 'x'.")
    print("This cross-section over a secondary domain 'u' is an Interval Type-2 (IT2) fuzzy set.")
    print("Its uncertainty is captured by upper and lower membership functions. We will formulate one of these bounds using a Gaussian function.\n")
    
    # --- 2. Define Components ---
    # Here we define the components of the mathematical formula to be printed.
    # The bound is represented as mu_B(u|x)
    # The right side is a Gaussian function: exp(-0.5 * ((u - c_x) / sigma_x)^2)
    lhs = "μ_B(u|x)"
    op = "="
    func = "exp"
    open_paren = "("
    close_paren = ")"
    coeff = "-0.5"
    mult = "*"
    term_u = "u"
    term_c = "c_x"
    term_sigma = "σ_x"
    div = "/"
    power = "^"
    exponent = "2"

    # --- 3. Assemble and Print the Final Equation ---
    print("The final mathematical formulation is assembled as follows:")
    print("-" * 60)
    
    # Printing each part of the equation as requested
    print("The Bound (as a function of u given x): ", lhs)
    print("Equals: ", op)
    print("The exponential function: ", func)
    print("With the argument: ", end="")
    
    # Explicitly printing the numbers and symbols for clarity
    sys.stdout.write(open_paren)
    # Highlight the first number as requested
    sys.stdout.write(f"<{coeff}> ")
    sys.stdout.write(f"{mult} {open_paren}")
    sys.stdout.write(f"({term_u} - {term_c}) {div} {term_sigma}")
    sys.stdout.write(close_paren)
    sys.stdout.write(f" {power} ")
    # Highlight the second number as requested
    sys.stdout.write(f"<{exponent}>")
    sys.stdout.write(close_paren)
    print("\n")

    print("-" * 60)

    # --- 4. Explain the Components ---
    print("Where the components represent:")
    print(f"  {lhs:<10}: The membership value of the bound (upper or lower) for the vertical cross-section.")
    print(f"  {term_u:<10}: The secondary input variable, whose domain is described by the cross-section.")
    print(f"  {term_c:<10}: The center (mean) of the Gaussian function, which is dependent on the primary input 'x'.")
    print(f"  {term_sigma:<10}: The spread (standard deviation) of the Gaussian, which is also dependent on 'x'.")
    print(f"  <{coeff}>, <{exponent}>: The numeric constants that define the Gaussian shape.")

# Execute the function to generate the output
formulate_it3_vertical_slice()

# The final mathematical formulation is provided below as the answer.
# Using unicode characters for a cleaner look.
final_answer_formula = "μ_B(u|x) = exp(-0.5 * ((u - c_x) / σ_x)^2)"

# The prompt requires enclosing the final answer in <<< >>>
# For example: <<<C>>>, <<<9.8>>>.
# In this case, the answer is the formula itself.
# To avoid printing this from the script, it's included as the final deliverable.
