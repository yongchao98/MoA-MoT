def identify_knockout_target_for_coumarate_production():
    """
    This script identifies the gene to knock out in Corynebacterium glutamicum
    to prevent the degradation of p-coumaric acid.
    """
    organism = "Corynebacterium glutamicum"
    product = "p-coumaric acid"
    intermediate = "protocatechuate (PCA)"
    enzyme = "protocatechuate 3,4-dioxygenase"
    target_gene = "pcaG"
    target_operon = "pcaGH"

    print(f"Goal: Prevent degradation of {product} in {organism}.")
    print("-" * 60)

    print("Step 1: The Degradation Funnel")
    print(f"In {organism}, aromatic compounds like {product} are broken down by being converted")
    print(f"into a central chemical intermediate, which is {intermediate}.")
    print("-" * 60)

    print("Step 2: The Point of No Return")
    print(f"This intermediate, {intermediate}, is then fed into the β-ketoadipate pathway.")
    print(f"The first, and most critical, step of this pathway is the cleavage of the aromatic")
    print(f"ring of {intermediate}. This is catalyzed by the enzyme {enzyme}.")
    print("-" * 60)

    print("Step 3: The Genetic Target")
    print(f"The enzyme {enzyme} is encoded by the '{target_gene}' and 'pcaH' genes.")
    print(f"These genes ({target_operon}) are the primary gateway for degrading {intermediate}.")
    print("-" * 60)

    print("Conclusion: The Knockout Strategy")
    print(f"To stop the degradation of {product}, we must block this gateway.")
    print(f"Knocking out the '{target_gene}' gene will disable the {enzyme}, creating a metabolic")
    print(f"block. This prevents the cell from breaking down the product, forcing it to accumulate.")
    print("\nTherefore, the recommended gene to knock out is:")
    print(f"Target Gene: {target_gene}")

if __name__ == '__main__':
    identify_knockout_target_for_coumarate_production()