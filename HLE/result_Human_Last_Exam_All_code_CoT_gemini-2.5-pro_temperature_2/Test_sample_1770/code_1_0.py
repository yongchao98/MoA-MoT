def find_knockout_target():
    """
    Identifies the target gene for knockout to prevent p-coumaric acid
    degradation in Corynebacterium glutamicum.
    """
    organism = "Corynebacterium glutamicum"
    product = "p-coumaric acid"
    problem = f"The desired product, {product}, is degraded by the host organism, {organism}."

    # Scientific literature has identified the specific degradation pathway.
    pathway = "cmt (coumarate) catabolic pathway"
    first_enzyme = "p-coumarate hydratase"
    target_gene_name = "cmtB"
    locus_tag = "NCgl1993" # In C. glutamicum ATCC 13032

    # Print the step-by-step reasoning
    print("--- Engineering Strategy: Preventing Product Degradation ---")
    print(f"Organism: {organism}")
    print(f"Target Product: {product}")
    print(f"Problem: {problem}\n")

    print("--- Analysis ---")
    print(f"1. Native Degradation Pathway: In {organism}, {product} is primarily degraded via the '{pathway}'.")
    print(f"2. First Enzymatic Step: The pathway is initiated by the enzyme '{first_enzyme}'.")
    print(f"3. Gene Encoding the Enzyme: This enzyme is encoded by the gene '{target_gene_name}' (locus tag: {locus_tag}).\n")

    print("--- Recommendation ---")
    print(f"To prevent the degradation of {product}, you should knock out the gene responsible for the first catabolic step.")
    print(f"The recommended target gene for knockout is: {target_gene_name} ({locus_tag}).")

if __name__ == "__main__":
    find_knockout_target()