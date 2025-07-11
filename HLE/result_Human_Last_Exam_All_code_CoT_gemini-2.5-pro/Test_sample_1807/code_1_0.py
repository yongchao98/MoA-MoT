def calculate_fst(pop1_counts, pop2_counts):
    """
    Calculates Fst for a biallelic locus using the formula:
    Fst = (Ht - Hs) / Ht
    where:
    Ht = Expected heterozygosity in the total population
    Hs = Average expected heterozygosity within subpopulations
    """
    # pop_counts is a list [count_allele_A, count_allele_B]
    n1 = sum(pop1_counts)
    n2 = sum(pop2_counts)
    
    # Allele frequencies in each subpopulation
    p1 = pop1_counts[0] / n1
    q1 = 1 - p1
    p2 = pop2_counts[0] / n2
    q2 = 1 - p2
    
    # Hs: Average expected heterozygosity within subpopulations
    # Hs = (2*p1*q1 * n1 + 2*p2*q2 * n2) / (n1 + n2) -> Simplified for our case
    hs1 = 2 * p1 * q1
    hs2 = 2 * p2 * q2
    Hs = (hs1 + hs2) / 2 # Assuming equal population sizes for simplicity
    
    # Ht: Expected heterozygosity in the total population
    total_n = n1 + n2
    p_total = (pop1_counts[0] + pop2_counts[0]) / total_n
    q_total = 1 - p_total
    Ht = 2 * p_total * q_total
    
    # Fst
    if Ht == 0:
        return 0 # No variation in the total population
    
    Fst = (Ht - Hs) / Ht
    return Fst

# --- Scenario: Marker on a Y-chromosome ---
# Population 1: Males (e.g., 50 individuals)
# They all have the Y-chromosome, so they have the marker.
# Let's call the 'presence' of the marker allele 'A' and 'absence' allele 'B'.
# Males are hemizygous, so we count chromosomes. 50 males have 50 Y chromosomes.
male_counts = [50, 0] # 50 'A' alleles, 0 'B' alleles

# Population 2: Females (e.g., 50 individuals)
# They are XX, they do not have a Y-chromosome, so they don't have the marker.
female_counts = [0, 50] # 0 'A' alleles, 50 'B' alleles

# Calculate Fst
fst_value = calculate_fst(male_counts, female_counts)

# Let's calculate it step-by-step for clarity in the output
p_male = male_counts[0] / sum(male_counts)
p_female = female_counts[0] / sum(female_counts)

Hs_male = 2 * p_male * (1 - p_male)
Hs_female = 2 * p_female * (1 - p_female)
Hs = (Hs_male + Hs_female) / 2

p_total = (male_counts[0] + female_counts[0]) / (sum(male_counts) + sum(female_counts))
Ht = 2 * p_total * (1 - p_total)

print("This script demonstrates why sex chromosomes cause high Fst between sexes.")
print("Consider a marker on the Y chromosome in a population of 50 males and 50 females.\n")

print("--- Allele Frequencies ---")
print(f"Frequency of marker in Males (p_male): {p_male:.2f}")
print(f"Frequency of marker in Females (p_female): {p_female:.2f}\n")

print("--- Heterozygosity Calculation (Fst = (Ht - Hs) / Ht) ---")
print(f"Step 1: Calculate average heterozygosity within subpopulations (Hs)")
print(f"   - Expected Heterozygosity in Males (Hs_male) = 2 * {p_male:.2f} * (1 - {p_male:.2f}) = {Hs_male:.2f}")
print(f"   - Expected Heterozygosity in Females (Hs_female) = 2 * {p_female:.2f} * (1 - {p_female:.2f}) = {Hs_female:.2f}")
print(f"   - Hs = ({Hs_male:.2f} + {Hs_female:.2f}) / 2 = {Hs:.2f}\n")

print(f"Step 2: Calculate heterozygosity in the total population (Ht)")
print(f"   - Total frequency of marker (p_total) = ({male_counts[0]} + {female_counts[0]}) / ({sum(male_counts)} + {sum(female_counts)}) = {p_total:.2f}")
print(f"   - Ht = 2 * {p_total:.2f} * (1 - {p_total:.2f}) = {Ht:.2f}\n")

print("Step 3: Calculate Fst")
print(f"   - Fst = (Ht - Hs) / Ht")
print(f"   - Fst = ({Ht:.2f} - {Hs:.2f}) / {Ht:.2f}")
print(f"   - Final Fst = {fst_value:.2f}\n")

print("An Fst of 1.0 indicates complete differentiation, which is expected for a Y-linked marker.")
print("This demonstrates that sex determination systems are a primary cause of genetic differentiation between sexes.")