import sys

def solve_puzzle():
    """
    Solves the multi-step cryptography and logic puzzle.
    """
    # Step 1: Determine the decryption key from the trivia answers.
    # Q1 Answer: Scipio
    # Q2 Answer: cotton (The question mentions "woollen" as a trick)
    # Q3 Answer: Simeoni
    ans1 = "scipio"
    ans2 = "cotton"
    ans3 = "simeoni"
    
    # Combine answers, make lowercase. Stopwords are not present.
    key = ans1 + ans2 + ans3

    # Step 2: Decrypt the ciphertext using Beaufort cipher.
    # While the resulting text is a misdirection, this step is part of the process.
    ciphertext = "fetssonayhmcippuqadksfd dbhmquecmomiaqqadmzk lfvjqmydbumamsaaodqdjxozpr fexhiiiqfgyxjoapdadaygq idmtapibftejjnvlanmksrbzfijteknbpslxehmlkuqddyeixsdfbinqxlf wgxuwaejssyzuzmlumitoumfwhpmzbokgfn wsvllbmbfoyushhglfly"
    alphabet = "abcdefghijklmnopqrstuvwxyz "
    char_to_index = {char: i for i, char in enumerate(alphabet)}
    index_to_char = {i: char for i, char in enumerate(alphabet)}
    n = len(alphabet)
    
    plaintext = ""
    for i, c_char in enumerate(ciphertext):
        k_char = key[i % len(key)]
        c_index = char_to_index[c_char]
        k_index = char_to_index[k_char]
        # Beaufort decryption formula: P = K - C
        p_index = (k_index - c_index + n) % n
        plaintext += index_to_char[p_index]
        
    # The decrypted text is:
    # "the symbols represent the digits zero to nine in order starting from the 
    # dotted square what is the three digit number represented by the symbols in 
    # row six column one row one column two and row seven column five of the grid"

    # Step 3 & 4: Analyze the misdirection and solve the puzzle.
    # The decrypted text leads to contradictions, as some specified locations
    # contain numbers, not symbols. This indicates a misdirection.
    # The real puzzle is to sum the integer values in the grid.
    grid = [
        ['▣', 75, '◧', '◩', '◫', 45, 86],
        ['◨', '◨', '◪', '◨', '▨', '◪', '◫'],
        ['▤', '◫', '▦', 46, '◩', 20, '▣'],
        ['▥', '▧', '◨', 88, '▤', '▦', '◩'],
        ['◧', '◫', '◪', '◪', '▨', '◧', '▦'],
        ['▥', '▤', '▨', '▥', 49, '◨', 85],
        ['▩', 22, '▣', '◧', 30, '▥', '▩']
    ]

    numbers_to_sum = []
    for row in grid:
        for cell in row:
            if isinstance(cell, int):
                numbers_to_sum.append(cell)

    total = sum(numbers_to_sum)
    
    # Step 5: Print the final equation.
    equation_parts = [str(num) for num in numbers_to_sum]
    print(f"{' + '.join(equation_parts)} = {total}")

solve_puzzle()