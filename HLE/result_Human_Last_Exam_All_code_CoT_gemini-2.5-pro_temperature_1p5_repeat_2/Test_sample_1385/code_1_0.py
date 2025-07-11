def solve_enclitic_order():
    """
    Determines and displays the correct order of Old Russian enclitics.
    
    The order is based on a well-established grammatical hierarchy:
    1. Conjunctions (e.g., 'бо')
    2. Emphatic/Interrogative Particles (e.g., 'же')
    3. Conditional Particles (e.g., 'бы')
    4. Pronominal Clitics (e.g., 'мя')
    5. Verbal Clitics (e.g., 'еси')
    """
    
    # The given enclitics are: бы, же, еси, бо, мя
    # We arrange them according to the linguistic rule.
    ordered_enclitics = ['бо', 'же', 'бы', 'мя', 'еси']
    
    # A placeholder for the first stressed word in a sentence.
    stressed_word = "слово" # "слово" means "word" in Russian
    
    # Combine the stressed word and the enclitics into a single list
    equation_parts = [stressed_word] + ordered_enclitics
    
    # Format and print the final equation showing the correct order
    print("The correct order is demonstrated by the following equation:")
    print(" + ".join(equation_parts))

solve_enclitic_order()