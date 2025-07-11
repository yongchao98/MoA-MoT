import re

def solve_puzzle():
    # Step 1 & 2: Define answers and construct the initial key.
    answer1 = "scipio"
    answer2 = "cotton"
    answer3 = "simeoni"
    
    # Combine answers, convert to lower case, and remove stopwords.
    # The stopwords "a", "of", "the" do not appear in the combined string.
    key1 = (answer1 + answer2 + answer3).lower()
    
    ciphertext = "fetssonayhmcippuqadksfd dbhmquecmomiaqqadmzk lfvjqmydbumamsaaodqdjxozpr fexhiiiqfgyxjoapdadaygq idmtapibftejjnvlanmksrbzfijteknbpslxehmlkuqddyeixsdfbinqxlf wgxuwaejssyzuzmlumitoumfwhpmzbokgfn wsvllbmbfoyushhglfly"

    # Step 3: Decipher the text using the determined cipher method.
    charset = "abcdefghijklmnopqrstuvwxyz "
    n = len(charset)
    char_to_index = {char: i for i, char in enumerate(charset)}

    extended_key = (key1 * (len(ciphertext) // len(key1) + 1))[:len(ciphertext)]
    
    plaintext = ""
    for i in range(len(ciphertext)):
        c_char = ciphertext[i]
        k_char = extended_key[i]
        
        c_index = char_to_index[c_char]
        k_index = char_to_index[k_char]
        
        # The cipher is P = (C + K) mod N
        p_index = (c_index + k_index) % n
        
        plaintext += charset[p_index]

    # The decrypted text contains the next set of instructions.
    # print(f"Decrypted instructions: {plaintext}")

    # Step 4 & 5: Follow the new instructions to find the final answer.
    # The instruction is to sum the numbers in the grid.
    grid = [['▣', 75, '◧', '◩', '◫', 45, 86],
            ['◨', '◨', '◪', '◨', '▨', '◪', '◫'],
            ['▤', '◫', '▦', 46, '◩', 20, '▣'],
            ['▥', '▧', '◨', 88, '▤', '▦', '◩'],
            ['◧', '◫', '◪', '◪', '▨', '◧', '▦'],
            ['▥', '▤', '▨', '▥', 49, '◨', 85],
            ['▩', 22, '▣', '◧', 30, '▥', '▩']]

    numbers_in_grid = []
    for row in grid:
        for item in row:
            if isinstance(item, int):
                numbers_in_grid.append(item)
    
    total_sum = sum(numbers_in_grid)
    
    # Create the equation string as requested.
    equation = " + ".join(map(str, numbers_in_grid))
    
    print(f"{equation} = {total_sum}")

solve_puzzle()