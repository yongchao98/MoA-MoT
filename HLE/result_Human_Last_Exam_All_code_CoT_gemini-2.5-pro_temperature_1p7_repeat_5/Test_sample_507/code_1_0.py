def find_similar_words():
    """
    This function presents two languages from different Asian cultures with no direct contact
    that share surprisingly similar words for "mom", "dad", and "broom".
    """
    language1 = "Malay (Austronesian family, Southeast Asia)"
    language2 = "Telugu (Dravidian family, South India)"
    
    # Words in Language 1 (Malay)
    mom1 = "Emak"
    dad1 = "Bapa"
    broom1 = "Penyapu"

    # Words in Language 2 (Telugu)
    mom2 = "Amma"
    dad2 = "Baabu" # Also Naana, but Baabu is used for father/respected elder
    broom2 = "Cheepuru"

    print("The two languages are Malay and Telugu.")
    print("-" * 40)
    print(f"Language 1: {language1}")
    print(f"Language 2: {language2}")
    print("-" * 40)
    
    print("Comparison of words:\n")
    
    # Printing the "equation" for Mom
    print(f"Word: Mom")
    print(f"{language1.split(' ')[0]}: {mom1}")
    print(f"{language2.split(' ')[0]}: {mom2}")
    print(f"Similarity: {mom1} ≈ {mom2}\n")

    # Printing the "equation" for Dad
    print(f"Word: Dad")
    print(f"{language1.split(' ')[0]}: {dad1}")
    print(f"{language2.split(' ')[0]}: {dad2}")
    print(f"Similarity: {dad1} ≈ {dad2}\n")

    # Printing the "equation" for Broom
    print(f"Word: Broom")
    print(f"{language1.split(' ')[0]}: {broom1}")
    print(f"{language2.split(' ')[0]}: {broom2}")
    print(f"Similarity: {broom1} ≈ {broom2}")

find_similar_words()
