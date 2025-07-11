def solve_field_theory_inner_product():
    """
    Derives and prints the expression for the inner product (φ, Dφ)
    in finite-temperature scalar field theory.
    """

    # Using string variables to represent mathematical expressions
    phi = "φ"
    d_tau_sq = "∂_τ²"
    nabla_sq = "∇²"
    m_sq = "m²"
    d4x = "d⁴x_E"  # d⁴x_E = dτ d³x
    integral_sign = "∫"
    operator_D = f"[ -{d_tau_sq} - {nabla_sq} + {m_sq} ]"

    print("Derivation of the inner product (φ, Dφ):")
    print("-" * 40)

    # Step 1: The Euclidean Action
    print("\nStep 1: The Euclidean Action (S_E)")
    print("The partition function Z is defined by a functional integral with the Euclidean action S_E:")
    print(f"Z = ∫ D{phi} exp(-S_E[{phi}])\n")
    print("For a free, neutral scalar field, the action is:")
    print(f"S_E[{phi}] = {integral_sign} {d4x} [ (1/2)(∂_μ {phi})² + (1/2){m_sq}{phi}² ]")
    print("where (∂_μ φ)² is shorthand for (∂_τ φ)² + (∇φ)² in Euclidean space.\n")

    # Step 2: Integration by Parts
    print("Step 2: Rewriting the Action via Integration by Parts")
    print("We can rewrite the kinetic term, ∫ (∂_μ φ)², using integration by parts.")
    print("Assuming boundary terms vanish (due to periodicity in time and φ -> 0 at spatial infinity):")
    print(f"{integral_sign} {d4x} (∂_μ {phi})²  = -{integral_sign} {d4x} {phi}(∂_μ²){phi}")
    print(f"where ∂_μ² = {d_tau_sq} + {nabla_sq} is the 4D Euclidean Laplacian.\n")

    # Step 3: Assembling the action to find D
    print("Step 3: Identifying the Operator D")
    print("Substituting the result from Step 2 back into the action S_E:")
    print(f"S_E[{phi}] = {integral_sign} {d4x} [ (1/2)(- {phi}(∂_μ²){phi}) + (1/2){m_sq}{phi}² ]")
    print(f"S_E[{phi}] = (1/2) {integral_sign} {d4x} {phi} [ -∂_μ² + {m_sq} ] {phi}")
    print("The action is now in the form S_E = (1/2)(φ, Dφ), which allows us to identify the operator D.\n")
    print(f"The operator D is: D = -∂_μ² + m² = {operator_D}\n")

    # Step 4: The Final Inner Product
    print("Step 4: Final Expression for the Inner Product (φ, Dφ)")
    print("The inner product (f, g) in this context is defined by the integral ∫ f(x)g(x) d⁴x_E.")
    print("Therefore, the inner product (φ, Dφ) is given by the full integral expression:")

    final_equation = f"(φ, Dφ) = {integral_sign} {d4x} {phi}(x) {operator_D} {phi}(x)"
    print("\n" + "="*len(final_equation))
    print(final_equation)
    print("="*len(final_equation))
    
    # Printing each component as requested
    print("\nComponents of the final equation:")
    print(f"  - Inner product structure: (φ, Dφ)")
    print(f"  - Integration measure: {integral_sign} {d4x}")
    print(f"  - Field: {phi}(x)")
    print(f"  - Operator D acting on the field: {operator_D} {phi}(x)")
    print(f"    - Time derivative component: -{d_tau_sq}")
    print(f"    - Spatial derivative component: -{nabla_sq}")
    print(f"    - Mass component: +{m_sq}")


if __name__ == "__main__":
    solve_field_theory_inner_product()
