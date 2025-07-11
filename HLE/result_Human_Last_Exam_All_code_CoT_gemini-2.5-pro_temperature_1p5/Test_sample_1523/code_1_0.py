import math

def display_inner_product_expression():
    """
    This function prints the mathematical expression for the inner product
    (ϕ, D_ϕ) for a neutral scalar field in finite-temperature field theory.
    """
    
    # --- Introduction ---
    print("In finite-temperature field theory, the partition function Z for a neutral scalar field ϕ is given by the functional integral:")
    print("  Z = ∫ Dϕ exp(-S_E[ϕ])")
    print("where S_E is the Euclidean action.")
    print("\nThe action is commonly written in the quadratic form S_E = (1/2)(ϕ, D_ϕ).")
    print("The inner product (ϕ, D_ϕ) is therefore twice the action.")
    print("-" * 60)

    # --- The Final Equation ---
    print("\nThe inner product in position space is given by the following integral:")
    
    # Define the components of the mathematical expression
    integral_part = "∫₀^β dτ ∫ d³x"
    field = "ϕ(τ, x)"
    
    # Define the terms of the differential operator D
    # The coefficients are the numbers in the equation
    coeff_time_deriv = -1
    op_time_deriv = "(∂/∂τ)²"
    
    coeff_space_deriv = -1
    op_space_deriv = "∇²"
    
    coeff_mass_term = 1
    op_mass_term = "m²"
    
    # Construct and print the full equation string
    # We display the numbers (coefficients) 1 and -1 explicitly in the final equation's breakdown.
    final_equation = f"(ϕ, D_ϕ) = {integral_part} {field} * [ -{op_time_deriv} - {op_space_deriv} + {op_mass_term} ] {field}"
    
    print("\n" + final_equation + "\n")
    
    # --- Explanation of Terms ---
    print("Where the components of the final equation are:")
    print(f"  - `(ϕ, D_ϕ)`: The inner product.")
    print(f"  - `{integral_part}`: The integral over Euclidean spacetime.")
    print(f"      - `τ`: Imaginary time, integrated from 0 to β (where β = 1/T, and T is temperature).")
    print(f"      - `x`: The 3 spatial coordinates, integrated over all space.")
    print(f"  - `{field}`: The neutral scalar field.")
    print(f"  - `D = -{op_time_deriv} - {op_space_deriv} + {op_mass_term}`: This is the Euclidean Klein-Gordon operator. Let's output each component and its numerical coefficient:")
    print(f"      - Time derivative part: `{op_time_deriv}` with a coefficient of {coeff_time_deriv}.")
    print(f"      - Spatial derivative part: `{op_space_deriv}` (the Laplacian) with a coefficient of {coeff_space_deriv}.")
    print(f"      - Mass part: `{op_mass_term}` with a coefficient of {coeff_mass_term}.")

if __name__ == '__main__':
    display_inner_product_expression()
    # The final answer is the mathematical expression derived and explained above.
    # As the result is a formula, we represent it as a string.
    final_answer = "(ϕ, D_ϕ) = ∫₀^β dτ ∫ d³x ϕ(τ, x) * [ -(∂/∂τ)² - ∇² + m² ] ϕ(τ, x)"
    # print(f"\n\n<<< {final_answer} >>>") # Suppressing this from final output per convention, as the script itself provides the formatted answer.