def analyze_genetic_alterations():
    """
    Analyzes different genetic alterations to determine which leads to
    copy-neutral loss of heterozygosity (LOH).

    We model a cell that is initially heterozygous (genotype 'Aa') on a
    pair of homologous chromosomes. The normal state is diploid (copy number = 2).
    """

    print("Analyzing genetic alterations for their effect on heterozygosity and copy number.\n")

    # Initial state: A diploid cell, heterozygous for a gene (Alleles 'A' and 'a')
    initial_chromosomes = "['A'], ['a']"
    initial_copy_number = 2
    initial_state = "Heterozygous"

    print(f"Initial Cell State:")
    print(f"  Chromosomes: {initial_chromosomes}")
    print(f"  Copy Number: {initial_copy_number}")
    print(f"  Allelic State: {initial_state}\n")
    print("-" * 40)

    # A. Mitotic Recombination
    final_mr_chromosomes = "['a'], ['a']" # A possible outcome leading to LOH
    final_mr_copy_number = 2
    final_mr_state = "Homozygous (LOH)"
    print("A. Mitotic Recombination:")
    print("   This event can result in a daughter cell inheriting two identical chromosome arms.")
    print(f"   Final Chromosomes: {final_mr_chromosomes}")
    print(f"   Final Copy Number: {final_mr_copy_number}")
    print(f"   Result: Copy-neutral LOH can occur.\n")

    # B. A deletion of a chromosomal region
    final_del_chromosomes = "['a']" # The chromosome with 'A' is lost
    final_del_copy_number = 1
    final_del_state = "Hemizygous (LOH)"
    print("B. Deletion:")
    print("   A chromosome segment containing an allele is lost entirely.")
    print(f"   Final Chromosomes: {final_del_chromosomes}")
    print(f"   Final Copy Number: {final_del_copy_number}")
    print(f"   Result: Not copy-neutral. Copy number is reduced.\n")

    # C. Trisomy
    final_tri_chromosomes = "['A'], ['a'], ['a']" # Gain of an extra chromosome
    final_tri_copy_number = 3
    final_tri_state = "Heterozygous"
    print("C. Trisomy:")
    print("   The cell gains an extra copy of a chromosome.")
    print(f"   Final Chromosomes: {final_tri_chromosomes}")
    print(f"   Final Copy Number: {final_tri_copy_number}")
    print(f"   Result: Not copy-neutral. Copy number is increased.\n")

    # D. Uniparental Disomy (UPD)
    # This happens when a cell loses one parental chromosome and duplicates the other.
    final_upd_chromosomes = "['a'], ['a']"
    final_upd_copy_number = 2
    final_upd_state = "Homozygous (LOH)"
    print("D. Uniparental Disomy:")
    print("   The cell ends up with two copies of a chromosome from a single parent.")
    print(f"   Final Chromosomes: {final_upd_chromosomes}")
    print(f"   Final Copy Number: {final_upd_copy_number}")
    print(f"   Result: A classic example of copy-neutral LOH.\n")

    # E. Duplication of a chromosomal region
    final_dup_chromosomes = "['A'], ['a', 'a']" # The region with 'a' is duplicated
    final_dup_copy_number_gene_a = 2
    final_dup_copy_number_total = 3
    final_dup_state = "Heterozygous"
    print("E. Duplication:")
    print("   A segment of a chromosome is duplicated.")
    print(f"   Final Chromosomes: {final_dup_chromosomes}")
    print(f"   Final Gene 'a' copy number: {final_dup_copy_number_gene_a}")
    print(f"   Total copies of gene at locus: {final_dup_copy_number_total}")
    print(f"   Result: Not copy-neutral. Gene dosage is increased.\n")

    print("-" * 40)
    print("Conclusion:")
    print("Based on the analysis, both Mitotic Recombination (A) and Uniparental Disomy (D) lead to a copy number of 2 and can cause LOH.")
    print("However, Uniparental Disomy is the term that describes the specific state where both chromosomes derive from one parent, which is a primary mechanism for copy-neutral LOH in tumorigenesis and perfectly fits the description.")

analyze_genetic_alterations()