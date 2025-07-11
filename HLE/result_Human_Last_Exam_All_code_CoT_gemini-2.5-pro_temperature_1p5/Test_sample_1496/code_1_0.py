import unicodedata

def solve_puzzle():
    """
    Solves the multi-step cryptography and logic puzzle.
    """
    # Step 1 & 2: Form the decryption key from the answers to the three questions.
    # Q1: In act II of Manfred by Lord Byron, which historical figure does the Chamois hunter mention to convey their refusal of fame?
    # A1: Scipio
    # Q2: In Chapter 30 of Nicholas Nickleby by Charles Dickens, what is the unusual woollen item Smike wears while performing as a tailor?
    # A2: Cotton (referring to the single white cotton stocking)
    # Q3: In The Tartar Steppe by Dino Buzzati, what is the surname of the officer who holds the title of Lieutenant-Colonel in Chapter 25 of the novel?
    # A3: Simeoni
    
    # Combine answers, make lowercase, remove stopwords (none present).
    key = "scipiocottonsimeoni"

    # Step 3: Decipher the text.
    # The cipher is a Vigenère cipher (P = C - K), despite the prompt naming Beaufort.
    ciphertext = "fetssonayhmcippuqadksfd dbhmquecmomiaqqadmzk lfvjqmydbumamsaaodqdjxozpr fexhiiiqfgyxjoapdadaygq idmtapibftejjnvlanmksrbzfijteknbpslxehmlkuqddyeixsdfbinqxlf wgxuwaejssyzuzmlumitoumfwhpmzbokgfn wsvllbmbfoyushhglfly"
    alphabet = "abcdefghijklmnopqrstuvwxyz "
    m = len(alphabet)
    char_to_int = {char: i for i, char in enumerate(alphabet)}
    int_to_char = {i: char for i, char in enumerate(alphabet)}

    # Repeat key to match ciphertext length
    full_key = (key * (len(ciphertext) // len(key) + 1))[:len(ciphertext)]

    plaintext = ""
    for i in range(len(ciphertext)):
        c_val = char_to_int[ciphertext[i]]
        k_val = char_to_int[full_key[i]]
        p_val = (c_val - k_val + m) % m
        plaintext += int_to_char[p_val]

    print("Decrypted Question:")
    print(plaintext)
    print("-" * 20)

    # Step 4, 5, 6: Solve the decrypted question by summing the numbers in the grid.
    grid = [['▣', 75, '◧', '◩', '◫', 45, 86],
            ['◨', '◨', '◪', '◨', '▨', '◪', '◫'],
            ['▤', '◫', '▦', 46, '◩', 20, '▣'],
            ['▥', '▧', '◨', 88, '▤', '▦', '◩'],
            ['◧', '◫', '◪', '◪', '▨', '◧', '▦'],
            ['▥', '▤', '▨', '▥', 49, '◨', 85],
            ['▩', 22, '▣', '◧', 30, '▥', '▩']]
    
    numbers_to_sum = []
    for row in grid:
        for item in row:
            if isinstance(item, int):
                numbers_to_sum.append(item)

    total_sum = sum(numbers_to_sum)
    
    equation_str = " + ".join(map(str, numbers_to_sum))
    
    print("Final Calculation:")
    print(f"{equation_str} = {total_sum}")


solve_puzzle()
<<<B>>>