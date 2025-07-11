def solve_heritability_problem():
    """
    This script logically derives the maximum variance a Polygenic Score (PGS)
    can explain based on the given broad-sense heritability.
    """

    # --- Given Information ---
    # The broad-sense heritability (H^2) is given.
    # H^2 = V(G) / V(P), where V(G) is total genetic variance and
    # V(P) is total phenotypic variance.
    H_squared = 0.5

    # --- Core Concepts & Definitions ---
    # 1. Total Genetic Variance (V(G)) is the sum of additive (V(A)),
    #    dominance (V(D)), and epistatic/interaction (V(I)) variance components.
    #    V(G) = V(A) + V(D) + V(I)
    #    Since variance components are non-negative, it is a mathematical certainty that:
    #    V(A) <= V(G)

    # 2. Narrow-sense heritability (h^2) is the proportion of variance from additive effects:
    #    h^2 = V(A) / V(P)

    # 3. A standard Polygenic Score (PGS) from GWAS primarily captures additive effects.
    #    Therefore, its maximum theoretical explanatory power (R_squared_PGS) is h^2.
    #    R_squared_PGS <= h^2

    # --- Step-by-Step Derivation ---
    print("Deriving the upper limit for the Polygenic Score's explanatory power:")
    print("="*70)

    print(f"1. We are given the broad-sense heritability, H^2 = {H_squared}")
    print("   This means V(G) / V(P) = 0.5, or that total genetic factors explain 50% of phenotypic variance.")
    print("-" * 70)

    print("2. The relationship between narrow-sense (h^2) and broad-sense (H^2) heritability is derived from their definitions.")
    print("   Since V(A) <= V(G), we can divide by V(P) to get:")
    print("   (V(A) / V(P)) <= (V(G) / V(P))")
    print("   This means: h^2 <= H^2")
    print("-" * 70)

    print("3. The maximum variance a standard PGS can explain is h^2. Let's call this R_squared_PGS.")
    print("   So, we have the following inequality chain:")
    print("   R_squared_PGS <= h^2 <= H^2")
    print("-" * 70)

    print("4. Now, we substitute the known value of H^2 into the final equation.")
    print(f"   Final Equation: R_squared_PGS <= h^2 <= {H_squared}")
    print("-" * 70)

    print("\nConclusion:")
    print(f"The analysis shows that the variance explained by the polygenic score (R_squared_PGS) must be less than or equal to {H_squared}, or {H_squared * 100}%.")
    print("Therefore, the statement 'The polygenic score can not explain more than 50% of the variance in the phenotype' is necessarily true.")

# Run the logical demonstration
solve_heritability_problem()