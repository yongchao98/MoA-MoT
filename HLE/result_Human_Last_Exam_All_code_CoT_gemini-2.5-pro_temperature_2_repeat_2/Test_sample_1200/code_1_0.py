def simulate_genetic_alterations():
    """
    Simulates different genetic alterations to demonstrate which one leads to
    copy-neutral loss of heterozygosity (LOH).
    """
    # Initial state: A normal heterozygous diploid cell
    # Alleles are different ('A', 'a'), and the copy number is 2.
    initial_alleles = ['A', 'a']
    initial_copy_number = 2

    print(f"Initial State: Alleles = {initial_alleles}, Copy Number = {initial_copy_number}\n")
    print("--- Simulating Alterations ---\n")

    # A. Mitotic Recombination
    # This can result in a daughter cell that is homozygous for one of the alleles
    # but maintains the overall chromosome number.
    mr_alleles = ['a', 'a']
    mr_copy_number = 2
    is_loh_mr = len(set(mr_alleles)) == 1
    is_copy_neutral_mr = mr_copy_number == initial_copy_number
    print(f"A. Mitotic Recombination:")
    print(f"   Result -> Alleles: {mr_alleles}, Copy Number: {mr_copy_number}")
    print(f"   - Loss of Heterozygosity: {is_loh_mr}")
    print(f"   - Copy-Neutral: {is_copy_neutral_mr}\n")


    # B. Deletion of a chromosomal region
    # One of the alleles is lost completely.
    del_alleles = ['a']
    del_copy_number = 1
    is_loh_del = len(set(del_alleles)) == 1 # Technically LOH as only one allele type remains
    is_copy_neutral_del = del_copy_number == initial_copy_number
    print(f"B. Deletion:")
    print(f"   Result -> Alleles: {del_alleles}, Copy Number: {del_copy_number}")
    print(f"   - Loss of Heterozygosity: {is_loh_del}")
    print(f"   - Copy-Neutral: {is_copy_neutral_del}\n")


    # C. Trisomy
    # An extra copy of the chromosome is gained.
    tri_alleles = ['A', 'a', 'a'] # Example where 'a' chromosome was duplicated
    tri_copy_number = 3
    is_loh_tri = len(set(tri_alleles)) == 1
    is_copy_neutral_tri = tri_copy_number == initial_copy_number
    print(f"C. Trisomy:")
    print(f"   Result -> Alleles: {tri_alleles}, Copy Number: {tri_copy_number}")
    print(f"   - Loss of Heterozygosity: {is_loh_tri}")
    print(f"   - Copy-Neutral: {is_copy_neutral_tri}\n")


    # D. Uniparental Disomy (via chromosome loss and reduplication)
    # Step 1: Lose the chromosome with allele 'A'.
    # Step 2: Duplicate the remaining chromosome with allele 'a'.
    upd_alleles = ['a', 'a']
    upd_copy_number = 2
    is_loh_upd = len(set(upd_alleles)) == 1
    is_copy_neutral_upd = upd_copy_number == initial_copy_number
    print(f"D. Uniparental Disomy:")
    print(f"   Result -> Alleles: {upd_alleles}, Copy Number: {upd_copy_number}")
    print(f"   - Loss of Heterozygosity: {is_loh_upd}")
    print(f"   - Copy-Neutral: {is_copy_neutral_upd}\n")


    # E. Duplication of a chromosomal region
    # A segment containing an allele is duplicated.
    dup_alleles = ['A', 'a', 'a'] # Example where the 'a' allele region is duplicated
    dup_copy_number = 3
    is_loh_dup = len(set(dup_alleles)) == 1
    is_copy_neutral_dup = dup_copy_number == initial_copy_number
    print(f"E. Duplication:")
    print(f"   Result -> Alleles: {dup_alleles}, Copy Number: {dup_copy_number}")
    print(f"   - Loss of Heterozygosity: {is_loh_dup}")
    print(f"   - Copy-Neutral: {is_copy_neutral_dup}\n")

    print("--- Conclusion ---")
    print("Both Mitotic Recombination (A) and Uniparental Disomy (D) result in a copy-neutral loss of heterozygosity.")
    print("However, Uniparental Disomy is a key mechanism and specific state that perfectly describes maintaining gene dosage (copy-neutral) despite the loss of one parental allele.")


if __name__ == '__main__':
    simulate_genetic_alterations()
