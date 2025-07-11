import textwrap

def solve_ecological_puzzle():
    """
    Provides an estimate for the age of males of two insect species at mating.
    """

    # Species 1: Tridactylophagus tartari
    species_1_name = "Tridactylophagus tartari"
    species_1_info = textwrap.dedent(f"""\
    For {species_1_name}:
    This insect is a member of the order Strepsiptera. The adult males of this order have a very short, ephemeral lifespan. They do not feed as adults and their sole purpose is to find a female and mate. Their total adult lifespan is typically less than six hours. Therefore, mating must occur very quickly after they emerge from their host. A best average estimate for the male's age at mating is a few hours.
    Estimated Age: 6 hours.""")

    # Species 2: Periplaneta americana
    species_2_name = "Periplaneta americana"
    species_2_info = textwrap.dedent(f"""\
    For {species_2_name}:
    This is the American cockroach. After its final molt to become an adult (eclosion), it needs a short period to become sexually mature. Research indicates that males can become sexually mature and capable of mating as early as two days after eclosion, even though they may have to wait for a female to be receptive. The question asks for the best estimate of the male's age, and the assumption is that they can mate "relatively shortly after eclosion".
    Estimated Age: 2 days.""")

    # Conclusion
    conclusion = textwrap.dedent("""\
    Comparing these findings with the answer choices, the pair "6 hours, two days" is the best fit.
    
    Final Answer Breakdown:
    - Estimated age for Tridactylophagus tartari male: 6 hours
    - Estimated age for Periplaneta americana male: 2 days
    """)

    print(species_1_info)
    print("\n" + "="*30 + "\n")
    print(species_2_info)
    print("\n" + "="*30 + "\n")
    print(conclusion)

solve_ecological_puzzle()