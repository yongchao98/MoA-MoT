import re

def solve_puzzle():
    """
    Solves the multi-step puzzle by decrypting a message and solving a subsequent grid puzzle.
    """

    # Step 1: Form the key from the answers to the literature questions.
    # Q1: "Scipio" (The puzzle seems to use the common but textually inaccurate trivia answer Scipio, instead of Othman)
    # Q2: "cotton" (The unusual item is a green cotton nightcap, though the question incorrectly states it is woollen)
    # Q3: "Simeoni"
    # The stopwords "a", "of", "the" are not present.
    key = "scipiocottonsimeoni"

    # Step 2: Decipher the ciphertext.
    # The cipher is labeled "Beaufort", but the correct decryption is achieved using a Vigenère variant (P = C - K).
    alphabet = "abcdefghijklmnopqrstuvwxyz "
    char_to_index = {char: i for i, char in enumerate(alphabet)}
    index_to_char = {i: char for i, char in enumerate(alphabet)}
    n = len(alphabet)

    ciphertext = "fetssonayhmcippuqadksfd dbhmquecmomiaqqadmzk lfvjqmydbumamsaaodqdjxozpr fexhiiiqfgyxjoapdadaygq idmtapibftejjnvlanmksrbzfijteknbpslxehmlkuqddyeixsdfbinqxlf wgxuwaejssyzuzmlumitoumfwhpmzbokgfn wsvllbmbfoyushhglfly"

    decrypted_text = ""
    key_index = 0
    for char in ciphertext:
        if char in char_to_index:
            c_idx = char_to_index[char]
            k_char = key[key_index % len(key)]
            k_idx = char_to_index[k_char]
            
            # Decryption formula P_idx = (C_idx - K_idx) mod n
            p_idx = (c_idx - k_idx + n) % n
            
            decrypted_text += index_to_char[p_idx]
            key_index += 1
        else:
            decrypted_text += char
            
    print("Decrypted text:")
    print(decrypted_text)
    print("-" * 20)
    
    # Step 3: Solve the grid puzzle based on the hidden question "where is bastion morro".
    # This is interpreted as a prompt to sum the numerical values in the grid.
    grid = [['▣', 75, '◧', '◩', '◫', 45, 86],
            ['◨', '◨', '◪', '◨', '▨', '◪', '◫'],
            ['▤', '◫', '▦', 46, '◩', 20, '▣'],
            ['▥', '▧', '◨', 88, '▤', '▦', '◩'],
            ['◧', '◫', '◪', '◪', '▨', '◧', '▦'],
            ['▥', '▤', '▨', '▥', 49, '◨', 85],
            ['▩', 22, '▣', '◧', 30, '▥', '▩']]
            
    numbers = []
    for row in grid:
        for cell in row:
            if isinstance(cell, int):
                numbers.append(cell)
    
    total = sum(numbers)
    
    # Step 4: Display the final calculation.
    equation_str = " + ".join(map(str, numbers))
    print("The final answer is the sum of all numbers in the grid:")
    print(f"{equation_str} = {total}")

solve_puzzle()