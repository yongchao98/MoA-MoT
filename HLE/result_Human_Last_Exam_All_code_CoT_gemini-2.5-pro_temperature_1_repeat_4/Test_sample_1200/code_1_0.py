def analyze_genetic_alterations():
    """
    Analyzes genetic alterations to find the one causing copy-neutral loss of heterozygosity.
    """
    # Define the criteria from the question
    target_criteria = {
        "copy_number": "neutral",
        "heterozygosity": "loss"
    }

    # Define the properties of each answer choice
    options = {
        "A": {
            "name": "Mitotic recombination",
            "description": "Exchange of genetic material during mitosis, which can lead to daughter cells that are homozygous for certain gene loci.",
            "copy_number": "neutral",
            "heterozygosity": "loss"
        },
        "B": {
            "name": "A deletion of a chromosomal region",
            "description": "A part of a chromosome is lost.",
            "copy_number": "loss",  # Results in one copy instead of two (monosomy)
            "heterozygosity": "loss"
        },
        "C": {
            "name": "Trisomy",
            "description": "Presence of an extra copy of a chromosome.",
            "copy_number": "gain", # Results in three copies instead of two
            "heterozygosity": "gain/maintained"
        },
        "D": {
            "name": "Uniparental disomy",
            "description": "An individual inherits two copies of a chromosome from one parent and no copy from the other.",
            "copy_number": "neutral", # Results in two copies, but from the same parent
            "heterozygosity": "loss"
        },
        "E": {
            "name": "Duplication of a chromosomal region",
            "description": "A segment of a chromosome is duplicated.",
            "copy_number": "gain", # Results in three or more copies of the region
            "heterozygosity": "gain/maintained"
        }
    }

    print("Analyzing which genetic alteration leads to copy-neutral loss of heterozygosity...\n")
    
    best_fit = None
    
    for key, properties in options.items():
        print(f"--- Evaluating Option {key}: {properties['name']} ---")
        
        is_copy_neutral = properties['copy_number'] == target_criteria['copy_number']
        causes_loh = properties['heterozygosity'] == target_criteria['heterozygosity']
        
        print(f"Is it copy-neutral? {'Yes' if is_copy_neutral else 'No'} (Actual: {properties['copy_number']})")
        print(f"Does it cause LOH? {'Yes' if causes_loh else 'No'} (Actual: {properties['heterozygosity']})")
        
        if is_copy_neutral and causes_loh:
            print("Verdict: This option matches the criteria.\n")
            if not best_fit: # Keep track of potential answers
                best_fit = []
            best_fit.append(key)
        else:
            print("Verdict: This option does not match the criteria.\n")
            
    print("--- Final Conclusion ---")
    print(f"Options that cause copy-neutral LOH: {best_fit}")
    print("Both Mitotic Recombination (A) and Uniparental Disomy (D) fit the criteria.")
    print("However, Uniparental Disomy (UPD) is the term that most accurately and fundamentally describes the state of having two homologous chromosomes from a single parent, which is the quintessential example of whole-chromosome copy-neutral LOH.")
    print("The process for UPD often involves losing a chromosome from one parent and duplicating the chromosome from the other parent, perfectly matching the description of 'maintaining the gene dosage despite allele deletion'.")
    print("\nTherefore, the best answer is D.")

analyze_genetic_alterations()
<<<D>>>