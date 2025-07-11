def solve_cohomology_question():
    """
    This function determines the minimal cohomology degree for significant
    non-trivial extensions and obstructions in semi-abelian categories.
    """

    # The degree for classifying general, non-split extensions.
    # H^1 classifies only the simpler 'split' extensions.
    # H^2 is the first group to classify all extensions.
    degree_for_nontrivial_extensions = 2

    # The degree for primary obstructions in classical obstruction theory.
    # When trying to lift a map, the first obstruction to its existence
    # lies in H^2.
    degree_for_primary_obstructions = 2

    # The question asks for the minimal degree at which both concepts
    # become significant. We take the minimum of the degrees where each
    # concept is first fundamentally addressed.
    minimal_degree = min(degree_for_nontrivial_extensions, degree_for_primary_obstructions)

    # Print the "equation" showing the components of our reasoning.
    print(f"The first degree to classify general non-trivial extensions is {degree_for_nontrivial_extensions}.")
    print(f"The first degree for primary obstructions is {degree_for_primary_obstructions}.")
    print(f"Therefore, the minimal significant degree is {minimal_degree}.")

solve_cohomology_question()