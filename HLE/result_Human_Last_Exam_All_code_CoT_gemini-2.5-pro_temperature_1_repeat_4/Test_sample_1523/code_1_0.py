import sys

def solve_inner_product():
    """
    This function explains and calculates the inner product (ϕ, D_ϕ)
    in finite-temperature field theory for a neutral scalar field.
    """

    # --- Introduction and Definitions ---
    print("Step 1: The Euclidean Action S[ϕ]")
    print("In finite-temperature field theory, the partition function is Z = ∫ Dϕ exp(-S[ϕ]).")
    print("For a free, massive, neutral scalar field ϕ, the Euclidean action S[ϕ] is:")
    action_initial = "S[ϕ] = ∫ d⁴x_E [ (1/2)(∂_μ ϕ)² + (1/2)m²ϕ² ]"
    print(action_initial)
    print("\nwhere x_E = (τ, x) is the Euclidean spacetime coordinate, and m is the mass.")
    print("-" * 40)

    # --- Identifying the Operator D ---
    print("\nStep 2: Identifying the Operator D")
    print("To find the operator D, we rewrite the action by integrating the kinetic term by parts.")
    print("∫ (∂_μ ϕ)² d⁴x_E = -∫ ϕ(∂_μ∂^μ)ϕ d⁴x_E  (discarding boundary terms)")
    print("This allows us to write the action in a quadratic form: S[ϕ] = (1/2)(ϕ, Dϕ).")
    action_operator_form = "S[ϕ] = (1/2) ∫ d⁴x_E ϕ [ -∂_μ∂^μ + m² ] ϕ"
    print(action_operator_form)
    print("\nFrom this, we can identify the operator D as the Euclidean Klein-Gordon operator:")
    operator_D = "D = -∂_μ∂^μ + m²"
    print(f"Operator D: {operator_D}")
    print("-" * 40)

    # --- Defining the Inner Product ---
    print("\nStep 3: The Inner Product (ϕ, Dϕ)")
    print("The inner product (f, g) is defined as ∫ f(x)g(x) d⁴x_E.")
    print("So, (ϕ, Dϕ) is:")
    inner_product_expression = "∫ d⁴x_E ϕ(x) [ -∂_μ∂^μ + m² ] ϕ(x)"
    print(f"(ϕ, Dϕ) = {inner_product_expression}")
    print("\nNote that this is exactly equal to twice the free-field action: (ϕ, Dϕ) = 2S[ϕ].")
    print("-" * 40)
    
    # --- Finite Temperature Considerations and Final Expression ---
    print("\nStep 4: Final Expression at Finite Temperature")
    print("At finite temperature T, the imaginary time τ is compactified on [0, β], where β = 1/T.")
    print("The integral ∫ d⁴x_E becomes ∫₀ᵝ dτ ∫ d³x.")
    print("The Laplacian ∂_μ∂^μ becomes ∂_τ² + ∇².")
    
    print("\nThe final equation for the inner product is constructed from the following components:")
    # Per the instruction, outputting each "number" (term) in the final equation
    print("  Integral: ∫₀ᵝ dτ ∫ d³x")
    print("  First field term: ϕ(τ, x)")
    print("  Operator term: [ -(∂_τ² + ∇²) + m² ]")
    print("  Second field term: ϕ(τ, x)")

    print("\nPutting all the components together, the final result is:")
    final_expression = "∫₀ᵝ dτ ∫ d³x ϕ(τ, x) [ -(∂_τ² + ∇²) + m² ] ϕ(τ, x)"
    print(final_expression)

solve_inner_product()