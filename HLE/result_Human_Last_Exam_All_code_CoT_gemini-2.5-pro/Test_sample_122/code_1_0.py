def solve_tapestry_mystery():
    """
    This function analyzes the metabolic pathways to determine the
    mutated enzyme and original color of an orange tapestry patch.

    The logic is based on deducing which combination of a single-enzyme
    mutation and an original color (yellow or blue) would result in
    a final mixture of red and yellow pigments, which appear orange.
    """

    # Final color is orange, which implies a mix of RED and YELLOW.

    # Source of RED: Yellow --(Enzyme A)--> Red.
    # This means Enzyme A must be functional and there must be a yellow substrate.

    # Source of YELLOW: Either original yellow pigment remains, or Blue --(Enzyme D)--> Yellow.

    # Let's test the hypothesis derived from the analysis:
    # Original Color: yellow
    # Mutated Enzyme: B
    # Process:
    # 1. Start with a yellow patch.
    # 2. Enzyme A (functional) converts some yellow to a red intermediate.
    # 3. Enzyme B (mutated) cannot convert the red intermediate further.
    # 4. The process stalls, leaving a mixture of the original yellow pigment
    #    and the newly created red intermediate.
    # 5. Yellow + Red = Orange. This matches the observation.

    # No other combination of a single mutation and original color produces orange.

    mutated_enzyme = "B"
    original_color = "yellow"

    # The required output format is <enzyme>-<colour>
    final_answer = f"{mutated_enzyme}-{original_color}"
    print(final_answer)

solve_tapestry_mystery()
<<<B-yellow>>>