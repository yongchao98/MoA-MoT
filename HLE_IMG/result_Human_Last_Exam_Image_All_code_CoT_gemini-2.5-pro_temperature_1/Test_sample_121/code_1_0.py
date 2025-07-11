import textwrap

def analyze_gene_evolution():
    """
    This function analyzes the provided bioinformatics evidence to determine the evolutionary relationship
    between genes in D. melanogaster and D. willistoni, and prints the most supported conclusion.
    """

    # The core task is to determine the relationship (ortholog/paralog) between the single
    # CG14879 gene in D. melanogaster and the two elements found in the corresponding
    # location in D. willistoni.

    # Reasoning Step 1: Evidence of Duplication and Separate Genes
    # D. melanogaster has one gene (CG14879) between 'pad' and 'Manf'.
    # D. willistoni has two elements in this same region. This implies a gene duplication event.
    # The text states RNA-Seq shows different expression patterns between the left and right elements.
    # This indicates they are two separate, independently regulated genes.
    reasoning1 = """
    The presence of two genomic elements in D. willistoni in a region where D. melanogaster has only one (CG14879) points to a gene duplication event in the D. willistoni lineage. The provided RNA-Seq data, which shows different expression levels and tissue-specific patterns between these two elements, strongly supports the conclusion that they are now two separate genes.
    """

    # Reasoning Step 2: Identifying the Ortholog via Parsimony
    # After a duplication, one gene is the ortholog (retains original function) and the other is a paralog.
    # The principle of parsimony suggests the gene most similar to the ancestral form is the ortholog.
    # D. melanogaster's CG14879 isoforms mostly have 2 Coding Sequences (CDSs).
    # The 'right element' in D. willistoni has predicted gene structures with 2 CDSs (e.g., 1a+2).
    # This structural conservation makes the 'right element' the most likely ortholog.
    reasoning2 = """
    To distinguish the ortholog (the gene retaining the original function) from the paralog (the copy), we use the principle of parsimony. The gene that most resembles the ancestral gene is the likely ortholog. The analysis shows that the 'right element' in D. willistoni can form a gene with two CDSs. This structure is the most parsimonious with the common isoforms of CG14879 in D. melanogaster, which also have two CDSs. Therefore, the 'right element' is the best candidate for the ortholog.
    """

    # Reasoning Step 3: Final Conclusion
    # The right element is the ortholog, the left is the paralog, and they are separate genes.
    # Answer choice 'P' captures all these points.
    reasoning3 = """
    Based on the evidence, the most logical conclusion is that the right element is the ortholog of CG14879, the left element is the resulting paralog from the duplication, and they function as separate genes. Answer choice 'P' perfectly summarizes this: 'The left element is most likely a paralog whereas the right element is most likely an ortholog of CG14879, because of the number of CDSs is more parsimonious between species for the right element. The elements are likely part of separate genes.'
    """

    print("Step-by-step reasoning:")
    print("-----------------------")
    print(textwrap.fill(reasoning1.strip(), width=80))
    print("\n" + textwrap.fill(reasoning2.strip(), width=80))
    print("\n" + textwrap.fill(reasoning3.strip(), width=80))
    print("\n-----------------------")
    print("The final answer is P.")

    # Final output format as requested
    print("<<<P>>>")

# Execute the analysis
analyze_gene_evolution()