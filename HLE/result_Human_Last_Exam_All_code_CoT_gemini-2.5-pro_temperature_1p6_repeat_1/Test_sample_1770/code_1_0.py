import textwrap

def find_knockout_target():
    """
    Identifies and explains the gene knockout target to prevent p-coumaric acid
    degradation in Corynebacterium glutamicum.
    """

    # Explanation of the metabolic issue
    explanation = """
    In Corynebacterium glutamicum, produced p-coumaric acid is degraded primarily through a catabolic pathway that starts with its activation to p-coumaroyl-CoA. This is the first committed step for degradation. By blocking this step, p-coumaric acid can no longer enter the pathway and will accumulate in the cell.
    """

    # The enzyme responsible
    enzyme = "4-coumarate-CoA ligase"
    ec_number = "6.2.1.12"
    
    # The gene encoding the enzyme in C. glutamicum ATCC 13032
    gene_locus_tag = "cg3125"
    gene_name = "phrC"

    # Print the findings
    print("To prevent the degradation of your product, p-coumaric acid, you should knock out the gene responsible for the first step in its catabolic pathway.")
    print("-" * 70)
    print(textwrap.fill(explanation.strip(), width=70))
    print("\n[+] The key enzyme is: {} (EC: {})".format(enzyme, ec_number))
    print("[+] This enzyme is encoded by the following gene in Corynebacterium glutamicum:")
    print("\n    Locus Tag: {}".format(gene_locus_tag))
    print("    Gene Name: {}".format(gene_name))
    
    print("-" * 70)
    print("\nConclusion: The gene that should be knocked out is {} ({}).".format(gene_locus_tag, gene_name))

if __name__ == '__main__':
    find_knockout_target()
