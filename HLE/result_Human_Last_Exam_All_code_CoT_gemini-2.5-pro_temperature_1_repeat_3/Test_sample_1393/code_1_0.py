def solve_puzzle():
    """
    This function solves the multi-layered puzzle by decoding Morse and Baudot codes,
    analyzing the resulting text, and determining the correct answer.
    """

    # --- Step 1: Decode the Morse Code Question ---
    morse_dict = {
        '.-': 'A', '-...': 'B', '-.-.': 'C', '-..': 'D', '.': 'E', '..-.': 'F',
        '--.': 'G', '....': 'H', '..': 'I', '.---': 'J', '-.-': 'K', '.-..': 'L',
        '--': 'M', '-.': 'N', '---': 'O', '.--.': 'P', '--.-': 'Q', '.-.': 'R',
        '...': 'S', '-': 'T', '..-': 'U', '...-': 'V', '.--': 'W', '-..-': 'X',
        '-.--': 'Y', '--..': 'Z', '-----': '0', '.----': '1', '..---': '2',
        '...--': '3', '....-': '4', '.....': '5', '-....': '6', '--...': '7',
        '---..': '8', '----.': '9', '.-.-.-': '.', '--..--': ',', '..--..': '?',
        '-..-.': '/', '-.--.': '(', '-.--.-': ')', '.-..-.': '"'
    }

    morse_question = "... . .-.. . -.-. - / - .... . / -.-. --- .-. .-. . -.-. - / --- .-. .. --. .. -. / --- ..-. / - .... . / ..-. --- .-.. .-.. --- .-- .. -. --. / ... . -. - . -. -.-. . / .-..-. - .... . / ..-. .- -. / -.-. .... .. . ..-. -.--.- - .. .- -. -.--.- - .- .. -. / .- -. -.. / -... .- -. -.. .. - ... / .- .-. . / .-- --- .-. - .... .-.. . ... ... / - --- / -- . -. - .. --- -. .-.-.- / .-- .. - .... / --- -. . / ... .. -. --. .-.. . / ... .-- --- .-. -.. --..-- / .. / -.-. .- -. / -... .-.. --- -.-. -.- / - .-. --- --- .--. ... / .- / -- .. .-.. .-.. .. --- -. .-.-.- .-..-. / - .... . / -.-. .... --- -.-. .. . / --. .. ...- . -. / -... . .-.. .-.. --- .-- / .. ... / . -. -.-. .-. -.-- .--. - . -.. / ..- ... .. -. --. / -... .- ..- -.. --- - / -.-. --- -.. ."
    
    decoded_question = "".join([
        "".join([morse_dict.get(c, '') for c in w.split(' ')]) + ' '
        for w in morse_question.split(' / ')
    ]).strip()

    print("Step 1: Decode the question from Morse Code.")
    print("---------------------------------------------")
    print(f"Decoded Question: {decoded_question}\n")

    # --- Step 2: Analyze the Decoded Question ---
    print("Step 2: Analyze the quote and find its source.")
    print("---------------------------------------------")
    print("The quote is: \"THE FAN CHIEF(TIAN)TAIN AND BANDITS ARE WORTHLESS TO MENTION. WITH ONE SINGLE SWORD, I CAN BLOCK TROOPS A MILLION.\"")
    print("This is a well-known boast by the general Guan Yu from the classic Chinese novel 'Romance of the Three Kingdoms'.")
    print("The question asks for the 'origin' of this sentence, and the answer choices are genres of Chinese Opera encrypted in Baudot code.\n")

    # --- Step 3: Decode the Baudot Code Answer Choices ---
    baudot_dict = {
        '00011': 'A', '11001': 'B', '01110': 'C', '01001': 'D', '00001': 'E',
        '01101': 'F', '11010': 'G', '10100': 'H', '00110': 'I', '01011': 'J',
        '01111': 'K', '10010': 'L', '11100': 'M', '01100': 'N', '11000': 'O',
        '10110': 'P', '10111': 'Q', '01010': 'R', '00101': 'S', '10000': 'T',
        '00111': 'U', '11110': 'V', '10011': 'W', '11101': 'X', '10101': 'Y',
        '10001': 'Z', '00100': ' '
    }
    
    choices = {
        'A': '01111 00111 01100 00100 10111 00111',
        'B': '10100 00111 00011 01100 11010 00100 11100 00001 00110 11101 00110',
        'C': '10101 00111 00001 00100 01011 00111',
        'D': '10101 00111 00100 01011 00111',
        'E': '01110 10100 00111 00011 01100 00100 01011 00111'
    }

    print("Step 3: Decode the Baudot code answer choices.")
    print("---------------------------------------------")
    decoded_choices = {}
    for choice, code in choices.items():
        decoded_text = "".join([baudot_dict.get(c, '?') for c in code.split(' ')])
        print(f"Choice {choice}: {code}  -->  {decoded_text}")
        decoded_choices[choice] = decoded_text
    
    # --- Step 4: Analyze Choices and Conclude ---
    print("\nStep 4: Analyze choices and select the answer.")
    print("---------------------------------------------")
    print("The decoded choices are different genres of Chinese Opera:")
    for choice, name in decoded_choices.items():
        print(f"  {choice}. {name}")
    print("\nThe story of Guan Yu is a staple in many of these opera forms. The question specifically asks for the 'origin'.")
    print("Kunqu Opera (KUN QU) is one of the oldest and most influential forms, dating back to the 16th century. It is often called the 'ancestor of a hundred operas,' including Beijing Opera.")
    print("Many classic stories, including plays about Guan Yu from 'Romance of the Three Kingdoms', have their theatrical roots in Kunqu.")
    print("Because of its age and foundational influence, Kunqu is the most appropriate answer for the 'origin' of the theatrical tradition of this story.")

solve_puzzle()
print("\n<<<A>>>")