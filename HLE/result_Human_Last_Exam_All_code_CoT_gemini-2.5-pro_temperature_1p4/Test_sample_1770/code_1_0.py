def suggest_knockout_target():
    """
    Identifies and explains the gene knockout target to prevent
    p-coumaric acid degradation in Corynebacterium glutamicum.
    """
    
    # Information about the metabolic engineering problem
    organism = "Corynebacterium glutamicum"
    product_of_interest = "p-coumaric acid"
    
    # Identified gene target for knockout
    gene_locus_tag = "cg2900"
    gene_common_name = "fic / couL"
    gene_function = "p-Coumarate:CoA ligase"
    
    # Explanation for the knockout strategy
    explanation = (
        "In {organism}, the primary pathway for degrading {product} begins with its "
        "activation to p-coumaroyl-CoA. This reaction is catalyzed by the enzyme "
        "{function}. The gene responsible for encoding this enzyme has the locus tag "
        "'{gene_locus_tag}' (also known as {gene_common_name}).\n\n"
        "To prevent the degradation of your desired product, you should knock out the "
        "'{gene_locus_tag}' gene. This will block the first step of the catabolic pathway, "
        "preventing the conversion of p-coumaric acid and leading to its accumulation."
    ).format(
        organism=organism,
        product=product_of_interest,
        function=gene_function,
        gene_locus_tag=gene_locus_tag,
        gene_common_name=gene_common_name
    )
    
    # Print the final recommendation
    print("--- Gene Knockout Recommendation ---")
    print(f"Organism: {organism}")
    print(f"Product: {product_of_interest}")
    print("-" * 34)
    print(explanation)
    print("\nFinal Answer: The gene to knock out is identified by the locus tag cg" + "2900")


if __name__ == "__main__":
    suggest_knockout_target()
