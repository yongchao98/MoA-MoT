import textwrap

def explain_bridges_experiment():
    """
    Explains the reasoning to determine the non-disjunction event.
    """
    
    explanation_title = "Step-by-step analysis of Bridges' Experiment Observation:"
    
    step1_title = "1. Offspring Genotype and Phenotype"
    step1_text = """
    The offspring is a male with genotype X0. This means he has one X chromosome and no Y chromosome.
    His phenotype (red eyes, miniature wings) is determined by the genes on this single X.
    Therefore, his specific X chromosome carries the alleles for red eyes (+) and miniature wings (m), making his genotype X+m0.
    """

    step2_title = "2. Tracing the Gametes"
    step2_text = """
    An X+m0 zygote is formed by the fusion of a paternal sperm and a maternal egg.
    Based on Bridges' classic findings, exceptional males like this inherit their X chromosome from the father.
    This means the parental gametes must have been:
    - Sperm from the father: X+m
    - Egg from the mother: 0 (a 'nullo-X' egg, containing no sex chromosome)
    """

    step3_title = "3. Identifying the Causal Event"
    step3_text = """
    A female with two X chromosomes (XX) normally produces eggs that all contain one X.
    The production of a '0' (nullo-X) egg is the result of an error during meiosis called non-disjunction.
    Non-disjunction is the failure of chromosomes to separate properly.
    """

    step4_title = "4. Pinpointing the Specific Event"
    step4_text = """
    This non-disjunction must have occurred in the female during egg formation (oogenesis).
    This rules out non-disjunction in the male (C), autosomal non-disjunction (D), and mutation (E).
    The event could happen in Meiosis I (failure of homologous chromosomes to separate) or Meiosis II (failure of sister chromatids to separate).
    The classic event described by Bridges, termed 'primary non-disjunction,' is the failure of the homologous pair of X chromosomes to separate. This occurs during Meiosis I.
    """

    conclusion_title = "5. Final Conclusion"
    conclusion_text = """
    The event that explains the production of the '0' egg from the female parent is the non-disjunction of her X chromosomes.
    Specifically, the canonical explanation for Bridges' results is non-disjunction during Meiosis I.
    This leads to the correct answer: Non-disjunction of the X chromosome in female meiosis I.
    """

    print(explanation_title)
    print("-" * len(explanation_title))
    
    print(f"\n{step1_title}")
    print(textwrap.dedent(step1_text).strip())

    print(f"\n{step2_title}")
    print(textwrap.dedent(step2_text).strip())

    print(f"\n{step3_title}")
    print(textwrap.dedent(step3_text).strip())
    
    print(f"\n{step4_title}")
    print(textwrap.dedent(step4_text).strip())

    print(f"\n{conclusion_title}")
    print(textwrap.dedent(conclusion_text).strip())

explain_bridges_experiment()