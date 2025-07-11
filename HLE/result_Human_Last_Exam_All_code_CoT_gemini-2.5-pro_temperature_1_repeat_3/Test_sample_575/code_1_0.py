import sys
import io

# Helper function to redirect print output
def execute_and_capture_output(func):
    old_stdout = sys.stdout
    sys.stdout = captured_output = io.StringIO()
    func()
    sys.stdout = old_stdout
    return captured_output.getvalue()

def solve():
    """
    Analyzes the relationship between heritability and polygenic scores to find the correct statement.
    """
    print("### Step-by-Step Analysis ###")
    print("\n---")
    print("Step 1: Defining the Key Concepts")
    print("---")
    print("Broad-Sense Heritability (H²): The proportion of the total variation in a phenotype that is due to all genetic factors combined.")
    print("  H² = (Additive Variance + Dominance Variance + Epistatic Variance) / Phenotypic Variance")
    print("\nNarrow-Sense Heritability (h²): The proportion of phenotypic variation due to only the additive effects of genes. This is the component of heritability that makes offspring resemble their parents.")
    print("  h² = Additive Variance / Phenotypic Variance")
    print("\nPolygenic Score (PGS): A score calculated for an individual based on their genotype at many genetic loci. Standard PGS models from GWAS are additive; they sum the estimated linear effects of SNPs to predict a trait.")
    
    print("\n---")
    print("Step 2: The Core Logical Relationships")
    print("---")
    print("From the definitions, we can establish two key points:")
    print("1. The maximum possible variance a standard, additive PGS can explain (let's call it R²_PGS) is equal to the narrow-sense heritability. This is because both are based on summing additive genetic effects.")
    print("   Therefore: R²_PGS ≤ h²")
    print("\n2. The additive genetic variance is just one component of the total genetic variance. Therefore, the additive variance must be less than or equal to the total genetic variance.")
    print("   This means: h² ≤ H²")

    print("\n---")
    print("Step 3: Applying the Given Information")
    print("---")
    print("The problem states that the broad-sense heritability is 0.5.")
    print("  H² = 0.5")
    print("\nWe can now combine all the information into a single chain of inequalities:")
    print("  R²_PGS ≤ h² ≤ H²")
    print("\nSubstituting the given value, we get our final equation:")
    print("  R²_PGS ≤ h² ≤ 0.5")

    print("\n---")
    print("Step 4: Evaluating the Answer Choices")
    print("---")
    print("Let's test each statement against our final equation (R²_PGS ≤ 0.5):")
    
    print("\nA. The polygenic score can not explain more than 50% of the variance in the phenotype.")
    print("   This statement says R²_PGS ≤ 0.5. Our logical derivation confirms this is ALWAYS true. This statement is necessarily true.")
    
    print("\nB. Given an arbitrarily large GWAS, the polygenic score will approach a variance explained of 50%.")
    print("   This would only be true if h² = 0.5. However, we only know that h² ≤ 0.5. If there are any non-additive effects (dominance or epistasis), then h² will be less than 0.5. So, this statement is NOT necessarily true.")

    print("\nC. Given an arbitrarily large GWAS, the polygenic score ... will not approach a variance explained of 50% due to gene-gene interactions and other non-linear effects...")
    print("   This assumes that non-additive effects MUST exist, meaning h² < 0.5. While plausible for many traits, it is not a logical necessity. In a theoretical case with no non-additive effects, h² would equal 0.5. So, this statement is NOT necessarily true.")

    print("\nD. The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5.")
    print("   Epigenetic effects are not formally part of the standard heritability model (Vg = Va + Vd + Vi) and are usually considered environmental. Their existence does not place a necessary mathematical constraint on h² relative to H². This is not necessarily true.")
    
    print("\nE. None of the other answer choices are correct.")
    print("   This is false because statement A is correct.")

    print("\n---")
    print("Conclusion")
    print("---")
    print("The only statement that logically must be true based on the provided information is A. The variance explained by a polygenic score is fundamentally limited by the heritability of the trait, and specifically for a standard additive PGS, it cannot exceed the broad-sense heritability.")

# Print the analysis
solve()
<<<A>>>