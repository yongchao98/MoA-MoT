import sys
import io

# Keep the original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# --- Problem Setup ---
# We define the variables based on the problem statement.
# Let's set Total Phenotypic Variance (Vp) to 100 for clear demonstration.
Vp = 100.0
# Broad-sense heritability (H2) is given as 0.5.
H2 = 0.5

# --- Calculations ---
# Total Genetic Variance (Vg) is the portion of Vp attributable to genetics.
Vg = H2 * Vp
# Environmental Variance (Ve) is the remaining portion of Vp.
Ve = Vp - Vg

# --- Explanation and Analysis ---
print("This script analyzes the problem by breaking down the concepts of heritability.")
print("-" * 60)
print("1. Defining the Core Components:")
print(f"Based on the problem, we have a Broad-sense Heritability (H_sq) of {H2}.")
print(f"Let's assume a Total Phenotypic Variance (Vp) of {Vp} for this example.")
print(f"The Total Genetic Variance (Vg) is Vp * H_sq = {Vp} * {H2} = {Vg}.")
print(f"The Environmental Variance (Ve) is Vp - Vg = {Vp} - {Vg} = {Ve}.")
print("-" * 60)

print("2. Understanding Genetic Variance and Polygenic Scores (PGS):")
print("Total Genetic Variance (Vg) is composed of:")
print("  - Additive Variance (Va)")
print("  - Dominance Variance (Vd)")
print("  - Epistatic (Interaction) Variance (Vi)")
print("So, the full equation for Vg is: Vg = Va + Vd + Vi.")
print("\nA standard Polygenic Score (PGS) from GWAS data is a linear model. It primarily captures the Additive Variance (Va).")
print("The maximum variance a perfect, linear PGS can explain is the Narrow-sense Heritability (h_sq), where h_sq = Va / Vp.")
print("-" * 60)

print("3. Evaluating Each Statement:\n")

# --- Analysis of Statement A ---
print(">>> Analyzing Statement A: The polygenic score can not explain more than 50% of the variance in the phenotype.")
print("A PGS is built from genetic data. Therefore, the variance it explains must be a part of the Total Genetic Variance (Vg).")
print("The maximum proportion of variance any genetic predictor can explain is the Broad-sense Heritability (H_sq).")
print(f"The calculation is: (Maximum Explainable Genetic Variance) / Vp <= Vg / Vp")
print(f"In this problem, Vg / Vp = {Vg} / {Vp} = {H2}.")
print("Therefore, the PGS can explain at most 50% of the phenotypic variance. It cannot explain the 50% of variance caused by the environment (Ve).")
print("CONCLUSION: Statement A is necessarily TRUE.\n")

# --- Analysis of Statements B and C ---
print(">>> Analyzing Statements B and C:")
print("Statement B claims the PGS will approach 50% variance explained.")
print("Statement C claims it will NOT approach 50% variance explained due to non-linear effects.")
print("Let's test two possible scenarios for the composition of Vg:\n")

# Scenario 1 for B/C: All genetic variance is additive
print("Scenario 1: Assume all genetic variance is purely additive.")
Va_1 = Vg
Vd_1 = 0.0
Vi_1 = 0.0
h2_1 = Va_1 / Vp
print(f"  Here, Vg ({Vg}) = Va ({Va_1}) + Vd ({Vd_1}) + Vi ({Vi_1}).")
print(f"  The variance explained by an ideal linear PGS (h_sq) would be: Va / Vp = {Va_1} / {Vp} = {h2_1}.")
print("  In this case, the PGS *does* approach 50%. So, Statement B is true and C is false.\n")

# Scenario 2 for B/C: Some genetic variance is non-additive
print("Scenario 2: Assume some genetic variance is non-additive.")
Va_2 = Vg * 0.7  # 70% of Vg is additive
Vd_2 = Vg * 0.2  # 20% is dominance
Vi_2 = Vg * 0.1  # 10% is epistatic
h2_2 = Va_2 / Vp
print(f"  Here, Vg ({Vg}) = Va ({Va_2}) + Vd ({Vd_2}) + Vi ({Vi_2}).")
print(f"  The variance explained by an ideal linear PGS (h_sq) would be: Va / Vp = {Va_2:.1f} / {Vp} = {h2_2}.")
print("  In this case, the PGS approaches 35%, which is NOT 50%. So, Statement B is false and C is true.\n")

print("CONCLUSION: Since the truth of B and C depends on the unknown composition of Vg, neither is *necessarily* true.\n")


# --- Analysis of Statement D ---
print(">>> Analyzing Statement D: The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5.")
print("Narrow-sense heritability (h_sq = Va/Vp) is already less than or equal to Broad-sense heritability (H_sq = Vg/Vp = 0.5).")
print("Epigenetic effects that are heritable can contribute to genetic variance. However, they don't *necessarily* contribute only to non-additive components (Vd or Vi). An epigenetic mark could have a simple additive effect.")
print("Therefore, we cannot be certain that h_sq would be strictly less than 0.5 just because epigenetic effects exist.")
print("CONCLUSION: Statement D is not necessarily true.")
print("-" * 60)

# Final output
sys.stdout = original_stdout
# Get the content from the buffer
output = captured_output.getvalue()
print(output)
print("FINAL CONCLUSION: Based on the analysis, only Statement A is a necessary truth derived from the fundamental definitions of heritability.")

<<<A>>>