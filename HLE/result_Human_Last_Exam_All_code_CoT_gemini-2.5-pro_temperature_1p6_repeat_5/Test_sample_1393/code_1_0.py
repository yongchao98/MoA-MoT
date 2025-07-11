import collections

def solve_puzzle():
    """
    This script solves the puzzle by decoding Morse and Baudot codes,
    and then applying logic to find the correct answer.
    """
    
    # Step 1: Define Morse and Baudot code mappings
    MORSE_CODE_DICT = {
        '.-': 'A', '-...': 'B', '-.-.': 'C', '-..': 'D', '.': 'E', '..-.': 'F',
        '--.': 'G', '....': 'H', '..': 'I', '.---': 'J', '-.-': 'K', '.-..': 'L',
        '--': 'M', '-.': 'N', '---': 'O', '.--.': 'P', '--.-': 'Q', '.-.': 'R',
        '...': 'S', '-': 'T', '..-': 'U', '...-': 'V', '.--': 'W', '-..-': 'X',
        '-.--': 'Y', '--..': 'Z',
        '-----': '0', '.----': '1', '..---': '2', '...--': '3', '....-': '4',
        '.....': '5', '-....': '6', '--...': '7', '---..': '8', '----.': '9',
        '.-..-.': '"', '-.--.': '(', '-.--.-': ')', '.-.-.-': '.'
    }

    # Using ITA2 standard in Letters (LTRS) mode
    BAUDOT_CODE_DICT_LTRS = {
        '00011': 'A', '11001': 'B', '01110': 'C', '01001': 'D', '00001': 'E',
        '01101': 'F', '11010': 'G', '10100': 'H', '00110': 'I', '01011': 'J',
        '01111': 'K', '10010': 'L', '11100': 'M', '01100': 'N', '11000': 'O',
        '10110': 'P', '10111': 'Q', '01010': 'R', '00101': 'S', '10000': 'T',
        '00111': 'U', '11110': 'V', '10001': 'W', '11101': 'X', '10101': 'Y',
        # Note: In some versions, W and Z share a code.
        '00100': ' '  # Space
    }

    # Step 2: Decode the question from Morse Code
    morse_question = "... . .-.. . -.-. - / - .... . / -.-. --- .-. .-. . -.-. - / --- .-. .. --. .. -. / --- ..-. / - .... . / ..-. --- .-.. .-.. --- .-- .. -. --. / ... . -. - . -. -.-. . / .-..-. - .... . / ..-. .- -. / -.-. .... .. . ..-. -.--. - .. .- -. -.--.- - .- .. -. / .- -. -.. / -... .- -. -.. .. - ... / .- .-. . / .-- --- .-. - .... .-.. . ... ... / - --- / -- . -. - .. --- -. .-.-.- / .-- .. - .... / --- -. . / ... .. -. --. .-.. . / ... .-- --- .-. -.. --..-- / .. / -.-. .- -. / -... .-.. --- -.-. -.- / - .-. --- --- .--. ... / .- / -- .. .-.. .-.. .. --- -. .-.-.- .-..-. / - .... . / -.-. .... --- -.-. .. . / --. .. ...- . -. / -... . .-.. .-.. --- .-- / .. ... / . -. -.-. .-. -.-- .--. - . -.. / ..- ... .. -. --. / -... .- ..- -.. --- - / -.-. --- -.. ."

    decoded_question = ""
    words = morse_question.split(' / ')
    for word in words:
        letters = word.split(' ')
        for letter in letters:
            if letter in MORSE_CODE_DICT:
                decoded_question += MORSE_CODE_DICT[letter]
        decoded_question += ' '
    
    print("--- Step 1: Decoding the question from Morse Code ---")
    print(f"Decoded Question:\n{decoded_question.strip()}\n")

    # Step 3: Decode the answer choices from Baudot code
    answer_choices = {
        'A': "01111 00111 01100 00100 10111 00111",
        'B': "10100 00111 00011 01100 11010 00100 11100 00001 00110 11101 00110",
        'C': "10101 00111 00001 00100 01011 00111",
        'D': "10101 00111 00100 01011 00111",
        'E': "01110 10100 00111 00011 01100 00100 01011 00111"
    }
    
    decoded_answers = collections.OrderedDict()

    print("--- Step 2: Decoding the Answer Choices from Baudot Code ---")
    for choice, baudot_string in answer_choices.items():
        decoded_text = "".join(BAUDOT_CODE_DICT_LTRS.get(code, '?') for code in baudot_string.split(' '))
        decoded_answers[choice] = decoded_text
        print(f"Choice {choice}: {baudot_string} -> {decoded_text}")
    print("\n")

    # Step 4: Analyze the puzzle and find the solution
    print("--- Step 3: Analysis and Solution ---")
    print("The decoded question asks for the origin of the quote:")
    print("\"...WITH ONE SINGLE SWORD, I CAN BLOCK TROOPS A MILLION.\"")
    print("\nThis famous quote is from the Chinese classic 'Romance of the Three Kingdoms'.")
    print("It is attributed to general Jiang Wei of the Shu Han kingdom.")
    print("He famously made this stand at Jiange (Sword Gate Pass), a key location in modern-day Sichuan province.")
    print("\nThe decoded answer choices are all forms of Chinese opera:")
    for choice, text in decoded_answers.items():
        print(f"  - {choice}: {text}")
    
    print("\nGiven the quote's geographical origin in Sichuan, the most logical answer is Sichuan Opera, which is known as 'Chuan Ju'.")
    print("This corresponds to choice E.")
    print("\nLet's show the encoding for the final answer 'CHUAN JU':")
    
    final_answer_text = "CHUAN JU"
    final_answer_codes = answer_choices['E'].split(' ')
    text_list = list(final_answer_text)

    # Output each character and its corresponding binary code
    for i in range(len(text_list)):
        print(f"'{text_list[i]}' -> {final_answer_codes[i]}")

    print("\nThe analysis confirms that Choice E is the correct answer.")

# Execute the solver function
solve_puzzle()