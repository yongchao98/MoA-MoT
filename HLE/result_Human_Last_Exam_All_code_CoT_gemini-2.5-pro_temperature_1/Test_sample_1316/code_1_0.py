def explain_bridges_experiment():
    """
    Explains the non-disjunction event in Bridges' Drosophila experiments
    that leads to an exceptional X0 male with red eyes and miniature wings.
    """

    # 1. Define the genetic context based on the problem.
    # The offspring's phenotype (red eyes, miniature wings) tells us the alleles on its X chromosome.
    # This X must have come from the father.
    father_X_chromosome = "X(w+, m)"  # w+ for red eyes, m for miniature wings
    father_genotype = (father_X_chromosome, "Y")

    # For a proper experimental cross, the mother would be homozygous for the alternative traits.
    mother_X_chromosome = "X(w, m+)" # w for white eyes, m+ for normal wings
    mother_genotype = (mother_X_chromosome, mother_X_chromosome)

    observed_offspring_genotype = "X(w+, m)0"
    observed_offspring_phenotype = "Red eyes, miniature wings, Male (sterile)"

    print("Analysis of Bridges' Non-Disjunction Experiment")
    print("-" * 60)
    print("This script deduces the cause of an exceptional offspring.")
    print(f"\nParental Cross Setup:")
    print(f"  - Mother Genotype: {mother_genotype[0]}/{mother_genotype[1]} (Phenotype: White eyes, Normal wings)")
    print(f"  - Father Genotype: {father_genotype[0]}/{father_genotype[1]} (Phenotype: Red eyes, Miniature wings)")

    print(f"\nObserved Exceptional Offspring:")
    print(f"  - Genotype: {observed_offspring_genotype}")
    print(f"  - Phenotype: {observed_offspring_phenotype}")
    print("-" * 60)

    # 2. Trace the origin of the gametes for the exceptional offspring.
    print("Tracing the Chromosomes for the Exceptional Male:")
    print(f"1. The male's phenotype is determined by his single X chromosome: {father_X_chromosome}.")
    print("2. This chromosome must be inherited from his father.")
    paternal_gamete = father_X_chromosome
    print(f"   => Paternal Gamete (Sperm): {paternal_gamete} (This is a normal gamete)")

    print("\n3. The male is X0, meaning he received NO sex chromosome from his mother.")
    maternal_gamete = "0"
    print(f"   => Maternal Gamete (Egg): {maternal_gamete} (This is an abnormal 'nullo-X' gamete)")
    print("-" * 60)

    # 3. Explain how the abnormal maternal gamete was formed.
    print("Explaining the Abnormal Egg Formation:")
    print("A '0' (nullo-X) egg is the result of non-disjunction in the mother during meiosis.")
    print("\nScenario: Non-disjunction in Female Meiosis I")
    print(f" - In a normal Meiosis I, the mother's homologous chromosomes ({mother_genotype[0]} and {mother_genotype[1]}) separate.")
    print(f" - In non-disjunction, they FAIL to separate and go to the same pole.")
    print(f" - This produces two types of eggs:")
    diplo_X_egg = f"{mother_genotype[0]}-{mother_genotype[1]}"
    nullo_X_egg = "0"
    print(f"   - One egg with both X's: '{diplo_X_egg}'")
    print(f"   - One egg with no X: '{nullo_X_egg}' <-- This is the one we need.")
    print("-" * 60)

    # 4. Show the final fertilization "equation".
    print("The Final Fertilization Event:")
    print("The abnormal maternal egg is fertilized by the normal paternal sperm.")
    print("\nEquation of Fertilization:")
    # The prompt asks to output each 'number' in the 'equation'.
    # We will output each component gamete and the resulting zygote.
    print(f"   Paternal Gamete: {paternal_gamete}")
    print(f"   Maternal Gamete: {nullo_X_egg}")
    print(f"          Result = {paternal_gamete} + {nullo_X_egg}")
    final_genotype = f"{paternal_gamete}0"
    print(f"   Final Genotype = {final_genotype}")

    print("\nThis matches the observed exceptional male. The error is non-disjunction of the X chromosome in female meiosis I.")

if __name__ == '__main__':
    explain_bridges_experiment()