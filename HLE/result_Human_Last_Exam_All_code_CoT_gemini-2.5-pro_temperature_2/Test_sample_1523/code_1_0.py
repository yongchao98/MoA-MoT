def solve_inner_product():
    """
    This function explains and prints the inner product (φ, D_φ)
    in finite-temperature scalar field theory.
    """

    # --- Define symbolic representations as strings for clarity ---
    phi = "φ"
    d_phi_d_tau = f"(∂{phi}/∂τ)"
    grad_phi = f"(∇{phi})"
    m = "m"

    # --- Introduction ---
    print("In finite-temperature field theory, the Euclidean action S[φ] for a neutral scalar field φ is:")
    print(f"S[φ] = (1/2) ∫ dτ d³x [ {d_phi_d_tau}² + {grad_phi}² + {m}²{phi}² ]")
    print("\nThis action can be expressed in the quadratic form S[φ] = (1/2)(φ, Dφ), where D is an operator.")
    print("By using integration by parts, the operator D is found to be:")
    print(f"D = -∂²/∂τ² - ∇² + {m}²")
    
    # --- The Inner Product Calculation ---
    print("\nThe inner product (φ, Dφ) is defined as ∫ d⁴x φ(Dφ).")
    print(f"(φ, Dφ) = ∫ dτ d³x {phi} [ (-∂²/∂τ² - ∇² + {m}²) {phi} ]")
    
    # --- Simplification and Final Result ---
    print("\nUsing integration by parts again (and dropping boundary terms), this simplifies.")
    print("The term with second derivatives becomes a term with first derivatives squared:")
    print(f"∫ dτ d³x [-{phi}(∂²/∂τ²){phi}] = ∫ dτ d³x {d_phi_d_tau}²")
    print(f"∫ dτ d³x [-{phi}(∇²){phi}] = ∫ dτ d³x {grad_phi}²")

    # Construct the final equation with explicit numerical coefficients as requested
    term1 = f"1 * {d_phi_d_tau}²"
    term2 = f"1 * {grad_phi}²"
    term3 = f"1 * {m}²{phi}²"
    final_equation = f"(φ, Dφ) = ∫ dτ d³x [ {term1} + {term2} + {term3} ]"

    print("\nTherefore, the final expression for the inner product is:")
    print(final_equation)
    print("\nNote: This is equal to 2 * S[φ].")

# Execute the function to print the solution
solve_inner_product()
