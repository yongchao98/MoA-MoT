def solve_heritability_problem():
    """
    This function explains the concepts of broad-sense and narrow-sense heritability
    to determine the correct answer for the given multiple-choice question.
    """

    print("### Step 1: Understanding Heritability Formulas ###")
    print("Phenotypic Variance (Vp) is the sum of Genetic Variance (Vg) and Environmental Variance (Ve).")
    print("Vp = Vg + Ve")
    print("Genetic Variance (Vg) can be broken down into Additive (Va), Dominance (Vd), and Epistatic (Vi) components.")
    print("Vg = Va + Vd + Vi\n")

    print("Broad-sense heritability (H^2) is the proportion of total variance due to ALL genetic factors.")
    print("H^2 = Vg / Vp = (Va + Vd + Vi) / (Va + Vd + Vi + Ve)\n")

    print("Narrow-sense heritability (h^2) is the proportion of total variance due to ONLY ADDITIVE genetic factors.")
    print("h^2 = Va / Vp = Va / (Va + Vd + Vi + Ve)\n")

    print("### Step 2: Analyzing the Rabbit Experiment ###")
    print("The problem states the rabbits' genetic variance is 'entirely additive'.")
    print("This means for the rabbits, Vd = 0 and Vi = 0.")
    print("Let's use the given H^2 of 0.75 and create a numerical example:")
    # Hypothetical values for the rabbit scenario
    Va_rabbit = 60
    Vd_rabbit = 0
    Vi_rabbit = 0
    Ve_rabbit = 20
    Vg_rabbit = Va_rabbit + Vd_rabbit + Vi_rabbit
    Vp_rabbit = Vg_rabbit + Ve_rabbit
    H2_rabbit = Vg_rabbit / Vp_rabbit
    h2_rabbit = Va_rabbit / Vp_rabbit

    print(f"Let Va = {Va_rabbit}, Vd = {Vd_rabbit}, Vi = {Vi_rabbit}, Ve = {Ve_rabbit}")
    print(f"Then Vg = {Va_rabbit} + {Vd_rabbit} + {Vi_rabbit} = {Vg_rabbit}")
    print(f"And Vp = {Vg_rabbit} + {Ve_rabbit} = {Vp_rabbit}")
    print(f"H^2 = Vg / Vp = {Vg_rabbit} / {Vp_rabbit} = {H2_rabbit}")
    print(f"h^2 = Va / Vp = {Va_rabbit} / {Vp_rabbit} = {h2_rabbit}")
    print("In this case, H^2 and h^2 are equal, as observed.\n")

    print("### Step 3: Finding the Cause for Differences between H^2 and h^2 ###")
    print("The difference between H^2 and h^2 arises when non-additive genetic variance (Vd or Vi) is greater than zero.")
    print("Let's model a scenario where H^2 and h^2 are different due to epistasis (Choice C).")
    # Hypothetical values for a different species with epistasis
    Va_species2 = 48
    Vd_species2 = 0
    Vi_species2 = 12  # Presence of epistatic variance
    Ve_species2 = 20
    Vg_species2 = Va_species2 + Vd_species2 + Vi_species2
    Vp_species2 = Vg_species2 + Ve_species2
    H2_species2 = Vg_species2 / Vp_species2
    h2_species2 = Va_species2 / Vp_species2

    print(f"Let Va = {Va_species2}, Vd = {Vd_species2}, Vi = {Vi_species2}, Ve = {Ve_species2}")
    print(f"Then Vg = {Va_species2} + {Vd_species2} + {Vi_species2} = {Vg_species2}")
    print(f"And Vp = {Vg_species2} + {Ve_species2} = {Vp_species2}")
    print(f"H^2 = Vg / Vp = {Vg_species2} / {Vp_species2} = {H2_species2:.2f}")
    print(f"h^2 = Va / Vp = {Va_species2} / {Vp_species2} = {h2_species2:.2f}")
    print("Here, the presence of epistatic variance (Vi > 0) causes H^2 to be greater than h^2.\n")

    print("### Step 4: Evaluating Answer Choices ###")
    print("A. Different environmental variances (E^2): Incorrect. Ve affects the denominator of both H^2 and h^2.")
    print("B. Misestimation of phenotypic variance: Incorrect. Vp is the denominator for both H^2 and h^2.")
    print("C. Presence of epistatic interactions: Correct. Epistasis (Vi) is a non-additive component included in H^2 but not in the numerator of h^2, causing them to differ.")
    print("D. Genetic linkage: Incorrect. Linkage is not the fundamental reason for the conceptual difference.")
    print("E. Variation in dominance genetic variance (D^2): While dominance (Vd) does cause a difference, the statement 'having no impact on h^2' is incorrect because Vd is part of the denominator (Vp). This makes choice C a better answer.\n")

    print("### Conclusion ###")
    print("The most accurate reason for a difference between broad-sense and narrow-sense heritability is the presence of non-additive genetic variance, such as epistasis or dominance. Choice C correctly identifies this.")

solve_heritability_problem()
print("<<<C>>>")