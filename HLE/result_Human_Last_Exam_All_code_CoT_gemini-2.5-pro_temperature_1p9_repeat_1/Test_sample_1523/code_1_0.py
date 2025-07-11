def display_inner_product_formula():
    """
    Calculates and displays the formula for the inner product (ϕ,D_ϕ)
    in the context of finite-temperature field theory for a neutral scalar field.
    """

    # --- Symbolic Representation of Equation Components ---
    # The inner product is defined by an integral over Euclidean spacetime.
    integral_symbol = "∫"
    spacetime_differential = "d⁴x"

    # The fields and operators involved.
    field = "ϕ(x)"
    operator_D_kinetic_part = "(-∂_μ∂^μ)"
    operator_D_mass_part = "m²"

    # --- Assembling the final equation string ---
    # The full operator D is the sum of the kinetic and mass parts.
    full_operator_D = f"[{operator_D_kinetic_part} + {operator_D_mass_part}]"

    # The inner product is (ϕ, D_ϕ) = ∫ d⁴x ϕ(x) * D * ϕ(x)
    equation = f"{integral_symbol} {spacetime_differential} {field} * {full_operator_D} * {field}"

    # --- Explanation of the physical context and derivation ---
    explanation = (
        "In finite-temperature quantum field theory, the partition function (Z) for a free, "
        "massive, neutral scalar field (ϕ) is given by a functional integral:\n"
        "Z = ∫ Dϕ exp(-S[ϕ])\n\n"
        "The Euclidean action, S[ϕ], for this field is:\n"
        "S[ϕ] = ∫ d⁴x * [ (1/2)(∂_μ ϕ)(∂^μ ϕ) + (1/2)m² ϕ² ]\n\n"
        "Using integration by parts on the kinetic term and assuming the boundary terms at infinity "
        "vanish, the action can be rewritten in a quadratic form:\n"
        "S[ϕ] = (1/2) * ∫ d⁴x ϕ(x) * [ -∂_μ∂^μ + m² ] * ϕ(x)\n\n"
        "This form is often expressed as S[ϕ] = (1/2) * (ϕ, D_ϕ), where D is the differential operator\n"
        "D = [ -∂_μ∂^μ + m² ], and the inner product (f, g) is defined as ∫ d⁴x f(x)g(x).\n\n"
        "Therefore, the inner product (ϕ, D_ϕ) is:"
    )

    # --- Print the full result to the user ---
    print(explanation)
    print("\n" + "="*60)
    print(f" (ϕ, D_ϕ) = {equation}")
    print("="*60 + "\n")

    # --- Detailed breakdown of symbols containing numbers, as requested ---
    print("Where the components of the equation are:")
    print(f"  -  {integral_symbol} {spacetime_differential}: This represents the integral over 4-dimensional Euclidean spacetime.")
    print("     The number '4' indicates the dimension of the spacetime (1 Euclidean time + 3 space).")
    print(f"  -  {field}: The neutral scalar field as a function of spacetime position x.")
    print(f"  -  {operator_D_kinetic_part}: The negative of the Euclidean d'Alembertian operator (Laplacian in 4D).")
    print(f"  -  {operator_D_mass_part}: The square of the field's mass, m.")
    print("     The number '2' is the exponent in the mass term.")


# Execute the function to provide the answer.
display_inner_product_formula()