def find_modifying_adjectives():
    """
    This function identifies and lists the adjectives modifying "mariti"
    in the given Latin sentence.
    
    Sentence: "Sed Maxentium suppositum ferunt arte muliebri tenere mariti animum 
               laborantis auspicio gratissimi partus coepti a puero."
    
    Analysis:
    - "mariti" is genitive singular ("of the husband").
    - "muliebri" is ablative singular, modifying "arte".
    - "laborantis" is genitive singular, modifying "mariti".
    - "gratissimi" is genitive singular, modifying "partus".
    - "coepti" is genitive singular, modifying "partus".
    
    Therefore, only one adjective modifies "mariti".
    """
    
    # The number of adjectives found to modify "mariti"
    count = 1
    
    # The list of the adjective(s)
    adjectives_list = ["laborantis"]
    
    print("The number of adjectives modifying 'mariti' is:")
    print(count)
    
    print("\nThe adjective is:")
    for adj in adjectives_list:
        print(adj)

find_modifying_adjectives()