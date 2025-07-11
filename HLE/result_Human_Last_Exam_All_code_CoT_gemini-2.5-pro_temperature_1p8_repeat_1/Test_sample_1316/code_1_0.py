def bridges_experiment_logic():
    """
    Explains the genetic cross leading to an exceptional X0 male in Drosophila.
    """
    # Define the parental genotypes that could produce the exceptional offspring.
    # The father must provide the X chromosome with the red-eye (+) and miniature-wing (m) alleles.
    father_genotype = "X+m Y"
    # The mother provides the nullo-gamete (O). Her own genotype isn't critical for this specific outcome,
    # but in a classic cross, she might be recessive for the father's traits (e.g., white-eyed, normal-winged).
    mother_phenotype = "Any (undergoes non-disjunction)"

    print("--- Tracing the Chromosomes for an Exceptional X0 Male ---")
    print(f"Observed Offspring Genotype: X0 (Male)")
    print(f"Observed Offspring Phenotype: Red eyes, miniature wings\n")

    # The phenotype tells us what the single X chromosome must be.
    offspring_x_chromosome = "X+m"
    print(f"This phenotype means the male's single X chromosome is: {offspring_x_chromosome}\n")

    # Determine the origin of the gametes based on Bridges' findings.
    print("Step 1: The male offspring inherits his single X chromosome from his father.")
    paternal_gamete = offspring_x_chromosome
    print(f"Father's contribution (sperm): {paternal_gamete}")

    print("\nStep 2: The offspring inherits no sex chromosome from his mother.")
    maternal_gamete = "O"
    print(f"Mother's contribution (egg): {maternal_gamete} (This is a 'nullo-gamete')")

    # Form the 'equation' for the zygote.
    zygote = paternal_gamete + maternal_gamete
    print("\n--- Final Equation of Fertilization ---")
    print(f"Paternal Gamete + Maternal Gamete = Zygote Genotype")
    # To satisfy the "output each number" constraint, we print the components.
    paternal_part = paternal_gamete
    maternal_part = maternal_gamete
    final_zygote = zygote
    print(f"        '{paternal_part}'   +        '{maternal_part}'        =    '{final_zygote}'")


    print("\nConclusion:")
    print("The mother produced an egg with no X chromosome ('O').")
    print("This is caused by a non-disjunction event of the X chromosomes during female meiosis.")
    print("Both non-disjunction in meiosis I or meiosis II in the female can lead to this outcome.")

bridges_experiment_logic()