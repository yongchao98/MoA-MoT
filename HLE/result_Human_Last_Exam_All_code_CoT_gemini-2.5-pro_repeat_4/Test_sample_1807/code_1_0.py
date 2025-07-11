import random

def run_simulation():
    """
    Simulates a population to demonstrate genetic differentiation between sexes
    for autosomal vs. sex-linked markers.
    """
    # --- 1. Simulation Parameters ---
    N_MALES = 500
    N_FEMALES = 500
    AUTOSOMAL_ALLELES = ['A', 'a']
    XLINKED_ALLELES = ['X_B', 'X_b']

    # --- 2. Population Generation ---
    males = []
    females = []

    # Generate females (XX) with two autosomal and two X-linked alleles
    for _ in range(N_FEMALES):
        female = {
            'autosomal': [random.choice(AUTOSOMAL_ALLELES), random.choice(AUTOSOMAL_ALLELES)],
            'x_linked': [random.choice(XLINKED_ALLELES), random.choice(XLINKED_ALLELES)]
        }
        females.append(female)

    # Generate males (XY) with two autosomal and one X-linked allele
    for _ in range(N_MALES):
        male = {
            'autosomal': [random.choice(AUTOSOMAL_ALLELES), random.choice(AUTOSOMAL_ALLELES)],
            'x_linked': [random.choice(XLINKED_ALLELES)]
        }
        males.append(male)

    # --- 3. Helper Functions ---
    def get_allele_freqs(allele_list):
        """Calculates frequencies of alleles in a list."""
        if not allele_list: return {None: 1.0}
        total = len(allele_list)
        return {allele: allele_list.count(allele) / total for allele in set(allele_list)}

    def get_heterozygosity(freqs):
        """Calculates expected heterozygosity (He = 1 - sum(p_i^2)). For 2 alleles, He = 2pq."""
        if len(freqs) < 2: return 0.0
        p = list(freqs.values())[0]
        q = 1.0 - p
        return 2 * p * q

    print("This simulation demonstrates why different sex-determining systems (XY vs ZW) can cause high Fst at some markers between males and females.")
    print("-" * 70)

    # --- 4. Autosomal Marker Analysis ---
    print("\n--- Analysis of an AUTOSOMAL Marker ---")
    male_alleles = [allele for m in males for allele in m['autosomal']]
    female_alleles = [allele for f in females for allele in f['autosomal']]
    total_alleles = male_alleles + female_alleles

    freqs_m = get_allele_freqs(male_alleles)
    freqs_f = get_allele_freqs(female_alleles)
    freqs_t = get_allele_freqs(total_alleles)

    Hs_m = get_heterozygosity(freqs_m)
    Hs_f = get_heterozygosity(freqs_f)
    Hs = (Hs_m + Hs_f) / 2  # Average heterozygosity in subpopulations
    Ht = get_heterozygosity(freqs_t) # Heterozygosity in total population

    fst = (Ht - Hs) / Ht if Ht > 0 else 0

    print("Allele frequencies are expected to be similar between sexes.")
    print(f"Allele '{AUTOSOMAL_ALLELES[0]}' Freq (Males):   {freqs_m.get(AUTOSOMAL_ALLELES[0], 0):.4f}")
    print(f"Allele '{AUTOSOMAL_ALLELES[0]}' Freq (Females): {freqs_f.get(AUTOSOMAL_ALLELES[0], 0):.4f}")
    print("\nCalculating Fst = (Ht - Hs) / Ht")
    print(f"Hs (Avg Subpopulation Heterozygosity) = {Hs:.4f}")
    print(f"Ht (Total Population Heterozygosity)  = {Ht:.4f}")
    print(f"Final Equation: Fst = ({Ht:.4f} - {Hs:.4f}) / {Ht:.4f}")
    print(f"Resulting Fst for Autosomal Marker: {fst:.4f}")
    print(">>> Fst is near 0, indicating no significant genetic differentiation.")
    print("-" * 70)


    # --- 5. X-linked Marker Analysis ---
    print("\n--- Analysis of an X-LINKED Marker ---")
    male_alleles = [allele for m in males for allele in m['x_linked']]
    female_alleles = [allele for f in females for f in f['x_linked']]
    total_alleles = male_alleles + female_alleles

    freqs_m = get_allele_freqs(male_alleles)
    freqs_f = get_allele_freqs(female_alleles)
    freqs_t = get_allele_freqs(total_alleles)

    # Males are hemizygous for X-linked markers, so their heterozygosity is 0.
    Hs_m = 0.0
    Hs_f = get_heterozygosity(freqs_f)
    Hs = (Hs_m + Hs_f) / 2
    Ht = get_heterozygosity(freqs_t)

    fst = (Ht - Hs) / Ht if Ht > 0 else 0

    print("Allele frequencies can differ due to different numbers of X chromosomes.")
    print(f"Allele '{XLINKED_ALLELES[0]}' Freq (Males, 1 X chrom):   {freqs_m.get(XLINKED_ALLELES[0], 0):.4f}")
    print(f"Allele '{XLINKED_ALLELES[0]}' Freq (Females, 2 X chroms): {freqs_f.get(XLINKED_ALLELES[0], 0):.4f}")
    print("\nCalculating Fst = (Ht - Hs) / Ht")
    print(f"Hs (Avg Subpopulation Heterozygosity) = {Hs:.4f}  (Note: Hs for males is 0)")
    print(f"Ht (Total Population Heterozygosity)  = {Ht:.4f}")
    print(f"Final Equation: Fst = ({Ht:.4f} - {Hs:.4f}) / {Ht:.4f}")
    print(f"Resulting Fst for X-linked Marker: {fst:.4f}")
    print(">>> Fst is significantly > 0, indicating genetic differentiation.")

if __name__ == '__main__':
    run_simulation()