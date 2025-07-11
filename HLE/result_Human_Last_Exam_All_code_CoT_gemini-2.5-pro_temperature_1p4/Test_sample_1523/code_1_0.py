def print_inner_product_expression():
    """
    Prints the expression for the inner product (ϕ, D_ϕ) for a neutral scalar field
    in the context of finite-temperature field theory.
    """

    # --- Using Unicode for clear mathematical notation ---
    phi = "\u03D5"
    D_op = "D"
    sum_symbol = "\u03A3"
    integral = "\u222B"
    partial = "\u2202"
    nabla = "\u2207"
    beta = "\u03B2"
    omega = "\u03C9"
    pi = "\u03C0"
    phi_tilde = f"{phi}\u0303"

    # The equation represents the quadratic part of the Euclidean action in momentum space.
    # Full action S_E = (1/2) * (ϕ, Dϕ)
    # The numbers in the equation are 2 (as in 2π) and d (for dimensions).
    # All powers are also 2.
    term1 = f"({phi}, {D_op}{phi})"
    sum_part = f"{sum_symbol}\u2099"  # Sum over n
    integral_part = f"{integral} d\u1d48p / ({2}{pi})\u1d48"  # ∫ dᵈp / (2π)ᵈ
    field_part = f"|{phi_tilde}(p, {omega}\u2099)|\u00b2"  # |ϕ̃(p, ωₙ)|²
    propagator_part = f"(p\u00b2 + {omega}\u2099\u00b2 + m\u00b2)"

    equation = f"{term1} = {sum_part} {integral_part} {field_part} {propagator_part}"

    # --- Print the formatted output ---
    print("In the functional integral representation of the partition function, the inner product (ϕ, D_ϕ) is the quadratic form in the exponent.")
    print("In momentum space, this is expressed as:")
    print("-" * 65)
    print(equation)
    print("-" * 65)
    print("\nWhere the terms are defined as:")
    print(f"  \u2022 {phi}: The neutral scalar field.")
    print(f"  \u2022 {D_op}: The Euclidean Klein-Gordon operator, where {D_op} = (-{partial}\u00b2_\u03c4 - {nabla}\u00b2 + m\u00b2).")
    print(f"  \u2022 {phi_tilde}(p, {omega}\u2099): The Fourier component of the field {phi}.")
    print(f"  \u2022 {sum_part}: A sum over all integers n (from -∞ to +∞).")
    print(f"  \u2022 p: The d-dimensional spatial momentum.")
    print(f"  \u2022 d: The number of spatial dimensions.")
    print(f"  \u2022 m: The mass of the field quantum.")
    print(f"  \u2022 {omega}\u2099: The bosonic Matsubara frequencies, defined as {omega}\u2099 = {2}{pi}n / {beta}.")
    print(f"  \u2022 {beta}: The inverse temperature, where {beta} = 1/T (in natural units where k_B = 1).")

if __name__ == "__main__":
    print_inner_product_expression()