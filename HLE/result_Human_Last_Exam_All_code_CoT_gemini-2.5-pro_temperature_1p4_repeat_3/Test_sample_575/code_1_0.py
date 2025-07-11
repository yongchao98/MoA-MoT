import numpy as np

def heritability_simulation():
    """
    This script simulates a population's genetics and phenotype to explore concepts of heritability.
    It demonstrates why a polygenic score's predictive power is capped by broad-sense heritability.
    """

    # --- Step 1: Define Parameters from the Problem ---
    # We are given a broad-sense heritability (H²) of 0.5. H² is the proportion of
    # phenotypic variance explained by ALL genetic factors (Vg / Vp).
    H2_broad_sense = 0.5

    # A standard Polygenic Score (PGS) primarily captures ADDITIVE genetic effects. The proportion
    # of variance explained by additive effects is the narrow-sense heritability (h² = Va / Vp).
    # By definition, h² must be less than or equal to H². We'll choose a value for h²
    # that is less than H² to demonstrate the role of non-additive genetic effects.
    h2_narrow_sense = 0.4

    # For simplicity in our simulation, let's assume the total variance of the phenotype (Vp) is 1.0.
    # From this, we can determine the variance of each component:
    total_phenotypic_variance = 1.0
    total_genetic_variance_Vg = H2_broad_sense * total_phenotypic_variance
    additive_genetic_variance_Va = h2_narrow_sense * total_phenotypic_variance
    non_additive_genetic_variance = total_genetic_variance_Vg - additive_genetic_variance_Va
    environmental_variance_Ve = total_phenotypic_variance - total_genetic_variance_Vg
    
    # Simulation parameters
    n_individuals = 20000
    n_snps = 500

    print("--- Analysis of the Statement ---")
    print(f"The total genetic contribution to phenotype variance (Broad-Sense Heritability, H²) is {H2_broad_sense*100}%.")
    print("A polygenic score (PGS) is a model built from genetic data.")
    print("Logically, a model based on genetics cannot explain more variance than genetics itself contributes.")
    print("Therefore, the maximum possible variance explained by any PGS is capped by H².")
    print(f"The equation that must hold is: Variance Explained by PGS <= H²")
    print(f"So, the PGS can explain at most {H2_broad_sense*100}% of the variance. Statement A is necessarily true.")
    print("The following simulation demonstrates this principle.")
    print("-" * 30)

    # --- Step 2: Simulate Genetic and Phenotypic Data ---
    
    # Simulate random genotypes (0, 1, or 2 for allele counts)
    allele_freqs = np.random.uniform(0.05, 0.5, n_snps)
    genotypes = np.random.binomial(2, allele_freqs, size=(n_individuals, n_snps))
    # Standardize genotypes to have mean=0 and variance=1 for stable calculations
    genotypes_std = (genotypes - np.mean(genotypes, axis=0)) / np.std(genotypes, axis=0)

    # Simulate true additive effects (betas) for each SNP
    true_additive_effects = np.random.normal(0, np.sqrt(additive_genetic_variance_Va / n_snps), n_snps)
    
    # Create the components for each individual
    additive_component = genotypes_std @ true_additive_effects
    non_additive_component = np.random.normal(0, np.sqrt(non_additive_genetic_variance), n_individuals)
    environmental_component = np.random.normal(0, np.sqrt(environmental_variance_Ve), n_individuals)
    
    # The final phenotype is the sum of all components
    phenotype = additive_component + non_additive_component + environmental_component
    
    # --- Step 3: Evaluate the 'Perfect' Polygenic Score ---
    # A perfect PGS from an infinitely large GWAS would perfectly identify the true additive component.
    pgs_prediction = additive_component
    
    # The variance explained by the PGS is the squared correlation between the phenotype and the PGS.
    variance_explained_by_pgs = np.corrcoef(phenotype, pgs_prediction)[0, 1]**2
    
    # For verification, calculate the broad-sense heritability in our simulated data.
    simulated_H2 = np.var(additive_component + non_additive_component) / np.var(phenotype)

    # --- Step 4: Display Results ---
    print("\n--- Simulation Results ---")
    print(f"Target Broad-Sense Heritability (H²): {H2_broad_sense:.3f}")
    print(f"Actual H² in simulated data: {simulated_H2:.3f}\n")
    print(f"Target Variance Explained by PGS (h²): {h2_narrow_sense:.3f}")
    print(f"Actual variance explained by 'perfect' PGS: {variance_explained_by_pgs:.3f}\n")
    print("--- Final Check ---")
    print("Is the variance explained by the PGS less than or equal to the broad-sense heritability?")
    final_equation = f"{variance_explained_by_pgs:.4f} <= {H2_broad_sense}"
    is_true = variance_explained_by_pgs <= H2_broad_sense
    print(f"The check is: {final_equation}  ... which is {is_true}.")
    print("The simulation confirms that the PGS explains a fraction of variance that does not exceed H².")

if __name__ == "__main__":
    heritability_simulation()
<<<A>>>