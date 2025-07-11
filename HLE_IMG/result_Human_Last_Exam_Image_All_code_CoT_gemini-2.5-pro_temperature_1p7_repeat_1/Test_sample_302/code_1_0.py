def get_seal_characters():
    """
    This function identifies and prints the simplified Chinese characters on the seal.
    """
    # The characters are identified from the seal, read top-to-bottom, right-to-left.
    # Traditional characters: 君, 子, 長, 壽
    # Simplified characters: 君, 子, 长, 寿
    
    char1_simplified = "君"
    char2_simplified = "子"
    char3_simplified = "长"
    char4_simplified = "寿"
    
    # The full phrase is "君子长寿" (jūn zǐ cháng shòu),
    # meaning "May the gentleman enjoy a long life."
    
    print(f"The character content on the seal in simplified Chinese is: {char1_simplified}{char2_simplified}{char3_simplified}{char4_simplified}")

get_seal_characters()