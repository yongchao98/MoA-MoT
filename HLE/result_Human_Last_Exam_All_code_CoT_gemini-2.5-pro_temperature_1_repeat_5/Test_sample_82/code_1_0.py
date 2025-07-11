def analyze_heritability():
    """
    Analyzes the statements about heritability and polygenic scores based on the problem description.
    """

    # --- Step 1: Define the parameters from the problem ---
    # For simplicity, we can set the total Phenotypic Variance (Vp) to 1.0.
    Vp = 1.0
    # The problem states broad-sense heritability (H^2) is 0.5.
    H2 = 0.5

    # Broad-sense heritability is the proportion of phenotypic variance explained by total genetic variance (Vg).
    # H^2 = Vg / Vp
    Vg = H2 * Vp

    # The remaining variance is due to the environment (Ve).
    # Vp = Vg + Ve
    Ve = Vp - Vg

    print("--- Step 1: Understanding the Basic Model ---")
    print(f"We assume Total Phenotypic Variance (Vp) = 1.0")
    print(f"Given Broad-Sense Heritability (H^2) = {H2}")
    print(f"Total Genetic Variance (Vg) is calculated as: Vg = H^2 * Vp = {H2} * {Vp} = {Vg}")
    print(f"Environmental Variance (Ve) is the remainder: Ve = Vp - Vg = {Vp} - {Vg} = {Ve}")
    print("This means all genetic factors combined explain 50% of the variance in the phenotype.\n")

    # --- Step 2: Decompose Genetic Variance ---
    # Total Genetic Variance (Vg) can be broken down into:
    # Va: Additive variance (captured by linear models like GWAS/PGS)
    # Vd: Dominance variance (non-linear effect at a single locus)
    # Vi: Epistatic variance (non-linear interactions between loci)
    # Vg = Va + Vd + Vi
    # Narrow-sense heritability (h^2) is defined as h^2 = Va / Vp.
    # A standard PGS from a large GWAS will explain variance approaching h^2.

    print("--- Step 2: Modeling a General Case for a Polygenic Trait ---")
    print("A 'polygenic trait' typically implies a complex genetic architecture, including non-additive effects.")
    print("Let's assume a plausible scenario where Vg (0.5) is split into components:")
    # This is an illustrative example where non-additive variance exists.
    Va_case = 0.35
    Vd_case = 0.10
    Vi_case = 0.05
    print(f"Example: Va={Va_case}, Vd={Vd_case}, Vi={Vi_case}. Note that {Va_case} + {Vd_case} + {Vi_case} = {Vg}")
    
    # In this case, the narrow-sense heritability (h^2) would be:
    h2_case = Va_case / Vp
    print(f"The predictable variance by a linear PGS (h^2) is: h^2 = Va / Vp = {Va_case} / {Vp} = {h2_case}")
    print(f"Note that h^2 ({h2_case}) is less than H^2 ({H2}). This gap is known as non-additive genetic variance.\n")


    # --- Step 3: Evaluate Each Statement ---

    print("--- Analysis of Statement A ---")
    print("Statement A: The polygenic score can not explain more than 50% of the variance in the phenotype.")
    print("A PGS is built from genetic data. The absolute maximum variance that can be explained by ANY model using genetics is the total genetic variance, Vg.")
    print(f"The maximum possible R^2 from genetics is Vg / Vp = {Vg} / {Vp} = {H2} (or 50%).")
    print("Therefore, a PGS can never explain MORE than 50%.")
    print("Conclusion: Statement A is NECESSARILY TRUE.\n")

    print("--- Analysis of Statement C ---")
    print("Statement C: ...the PGS... will not approach... 50% due to gene-gene interactions... and dominance.")
    print("A standard, linear PGS explains variance approaching h^2 = Va / Vp.")
    print("It will only approach 50% (H^2) if all genetic variance is additive (Va = Vg).")
    print("For any complex polygenic trait, it's expected that Vd or Vi > 0, which makes Va < Vg.")
    print(f"In our example, the PGS approaches {h2_case*100}%, not 50%, precisely because of Vd and Vi.")
    print("Conclusion: Statement C is considered TRUE in the standard context of polygenic traits.\n")

    print("--- Analysis of Statement B ---")
    print("Statement B: ...the PGS... will approach a variance explained of 50%.")
    print("This is the opposite of statement C. As explained, this would only happen in the specific case of purely additive genetics, which is not the general rule for polygenic traits.")
    print("Conclusion: Statement B is NOT necessarily true.\n")

    print("--- Analysis of Statement D ---")
    print("Statement D: The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5.")
    print("By definition, h^2 = Va/Vp is always less than or equal to H^2 = Vg/Vp, so h^2 <= 0.5.")
    print("This mathematical limit exists independently of epigenetics. Epigenetic effects that are not stably inherited are part of the environmental variance (Ve).")
    print("The existence of epigenetic variance does not force Va to be less than Vg.")
    print("Conclusion: Statement D is NOT necessarily true.\n")

    # --- Step 4: Final Conclusion ---
    print("--- Final Summary ---")
    print("Statement A is a hard ceiling based on the definition of broad-sense heritability and is always true.")
    print("Statement C accurately describes why the variance explained by a standard PGS is typically less than the broad-sense heritability for a complex trait.")
    print("Both A and C are correct statements under the standard interpretation of the terms.")

analyze_heritability()
<<<E>>>