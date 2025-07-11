def analyze_bridges_experiment():
    """
    Analyzes the formation of exceptional offspring in Bridges' Drosophila experiment
    to identify the specific non-disjunction event.
    """
    # Define parental genotypes and the target exceptional offspring
    female_parent = ('Xw', 'Xw')
    male_parent = ('X+', 'Y')
    target_male_offspring = 'X+O' # Genotype of the red-eyed, sterile male

    print(f"Investigating the origin of the exceptional male offspring with genotype: {target_male_offspring}")
    print("-" * 60)

    # Step 1: Determine the required gametes for the target offspring
    # The 'X+' must come from the father, as the mother only has 'Xw' alleles.
    required_sperm = 'X+'
    # To get an X0 genotype, the egg must have no sex chromosome.
    required_egg = 'O'

    print(f"Parental Cross: White-eyed Female ({''.join(female_parent)}) x Red-eyed Male ({''.join(male_parent)})")
    print(f"Formation of '{target_male_offspring}': Requires a '{required_sperm}' sperm and an '{required_egg}' egg.\n")

    # Step 2: Evaluate which meiotic error in the female could produce an 'O' egg.
    print("Evaluating potential non-disjunction events in the female:")

    # Scenario A: Non-disjunction in Meiosis I
    # The two homologous Xw chromosomes fail to separate. This produces both
    # n+1 (XwXw) and n-1 (O) gametes.
    female_gametes_meiosis1_ndj = ['XwXw', 'O']
    can_form_with_A = required_egg in female_gametes_meiosis1_ndj
    print(f"  A. Meiosis I Non-disjunction:")
    print(f"     - Generates gametes: {female_gametes_meiosis1_ndj}")
    print(f"     - Can it produce the required '{required_egg}' egg? {'Yes' if can_form_with_A else 'No'}")

    # Scenario B: Non-disjunction in Meiosis II
    # Failure of sister chromatids to separate. This can also produce n+1 and n-1 gametes.
    female_gametes_meiosis2_ndj = ['XwXw', 'O'] # Plus normal gametes from other meiotic events
    can_form_with_B = required_egg in female_gametes_meiosis2_ndj
    print(f"  B. Meiosis II Non-disjunction:")
    print(f"     - Generates gametes: {female_gametes_meiosis2_ndj} (among others)")
    print(f"     - Can it produce the required '{required_egg}' egg? {'Yes' if can_form_with_B else 'No'}\n")


    # Step 3: Use Bridges' complete findings to select the best explanation.
    # Bridges also found exceptional XwXwY white-eyed females.
    required_egg_for_exceptional_female = 'XwXw'
    print("Deciding between A and B using Bridges' complete observations:")
    print("Bridges also observed exceptional white-eyed females (XwXwY).")
    print(f"This requires an '{required_egg_for_exceptional_female}' egg from the mother.\n")
    print(f"The event of primary non-disjunction is defined as the failure of homologous chromosomes to separate.")
    print("This is, by definition, an error in Meiosis I.")
    print("This single mechanism (Meiosis I NDJ) explains the production of BOTH 'XwXw' and 'O' gametes, accounting for all of Bridges' exceptional offspring.")
    print("-" * 60)
    print("Conclusion: The most specific and complete answer is Non-disjunction in female meiosis I.")

analyze_bridges_experiment()