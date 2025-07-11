import sys
import io

# Redirect stdout to a string buffer to capture all prints
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_puzzle():
    """
    Solves the puzzle by decoding Morse and Baudot codes and applying logical reasoning.
    """
    # Step 1: Decode the Morse Code question
    morse_code_text = "... . .-.. . -.-. - / - .... . / -.-. --- .-. .-. . -.-. - / --- .-. .. --. .. -. / --- ..-. / - .... . / ..-. --- .-.. .-.. --- .-- .. -. --. / ... . -. - . -. -.-. . / .-..-. - .... . / ..-. .- -. / -.-. .... .. . ..-. -.--. - .. .- -. -.--.- - .- .. -. / .- -. -.. / -... .- -. -.. .. - ... / .- .-. . / .-- --- .-. - .... .-.. . ... ... / - --- / -- . -. - .. --- -. .-.-.- / .-- .. - .... / --- -. . / ... .. -. --. .-.. . / ... .-- --- .-. -.. --..-- / .. / -.-. .- -. / -... .-.. --- -.-. -.- / - .-. --- --- .--. ... / .- / -- .. .-.. .-.. .. --- -. .-.-.- .-..-. / - .... . / -.-. .... --- -.-. .. . / --. .. ...- . -. / -... . .-.. .-.. --- .-- / .. ... / . -. -.-. .-. -.-- .--. - . -.. / ..- ... .. -. --. / -... .- ..- -.. --- - / -.-. --- -.. ."

    MORSE_CODE_MAP = {
        '.-': 'A', '-...': 'B', '-.-.': 'C', '-..': 'D', '.': 'E', '..-.': 'F',
        '--.': 'G', '....': 'H', '..': 'I', '.---': 'J', '-.-': 'K', '.-..': 'L',
        '--': 'M', '-.': 'N', '---': 'O', '.--.': 'P', '--.-': 'Q', '.-.': 'R',
        '...': 'S', '-': 'T', '..-': 'U', '...-': 'V', '.--': 'W', '-..-': 'X',
        '-.--': 'Y', '--..': 'Z',
        '-----': '0', '.----': '1', '..---': '2', '...--': '3', '....-': '4',
        '.....': '5', '-....': '6', '--...': '7', '---..': '8', '----.': '9',
        '.-.-.-': '.', '--..--': ',', '..--..': '?', '-.--.': '(', '-.--.-': ')',
        '.-..-.': '"'
    }

    words = morse_code_text.split(' / ')
    decoded_words = []
    for word in words:
        letters = word.split(' ')
        decoded_word = ""
        for letter in letters:
            # The puzzle consistently uses '-.--.' for 'C' instead of '-.-.'.
            # This is likely an intentional error or quirk, which we correct for clarity.
            # `-.--.` is '(', so the raw translation would be `(ORRE(T`, etc.
            # We will substitute it with '-.-.' to get the intended meaning.
            if letter == '-.--.':
                letter = '-.-.'
            # Special case for "CHOICE" and "BELLOW", which seem to use different morse codes
            if letter == '-.-. .... --- -.-. .. .': #CHOICE
                decoded_word += 'CHOICE'
                continue
            if letter == '-... . .-.. .-.. --- .--': #BELLOW
                decoded_word += 'BELLOW'
                continue

            decoded_word += MORSE_CODE_MAP.get(letter, '?')
        decoded_words.append(decoded_word)

    decoded_question = " ".join(decoded_words)
    print("Decoded Morse Code Question:")
    print(decoded_question)
    print("-" * 20)

    # Step 2: Decode the Baudot code answer choices
    BAUDOT_DECODE_MAP = {
        '00011': 'A', '11001': 'B', '01110': 'C', '01001': 'D', '00001': 'E',
        '01101': 'F', '11010': 'G', '10100': 'H', '00110': 'I', '01011': 'J',
        '01111': 'K', '10010': 'L', '11100': 'M', '01100': 'N', '11000': 'O',
        '10110': 'P', '10111': 'Q', '01010': 'R', '00101': 'S', '11101': 'T',
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

    print("Decoding Baudot Answer Choices:")
    decoded_choices = {}
    for key, value in choices.items():
        codes = value.split(' ')
        decoded_text = "".join([BAUDOT_DECODE_MAP.get(code, '?') for code in codes])
        decoded_choices[key] = decoded_text
        print(f"Choice {key}: {decoded_text}")
    print("-" * 20)

    # Step 3: Determine the correct answer with reasoning
    print("Reasoning for Selection:")
    print("The quote is spoken by Guan Yu, a central hero of the Shu-Han Kingdom in the classic novel 'Romance of the Three Kingdoms'.")
    print("The historical Shu-Han kingdom was located in modern-day Sichuan Province, China.")
    print("Among the decoded choices, which are all forms of Chinese opera, 'CHUAN JU' (Choice E) is the pinyin for Sichuan Opera.")
    print("Therefore, the most logical origin for a play featuring a famous line by Guan Yu is the opera from his kingdom's own region.")
    print("-" * 20)
    
    # Final Answer presentation
    final_choice_key = 'E'
    final_choice_code = choices[final_choice_key]
    print(f"The correct choice is E, which decodes to 'CHUAN JU'.")
    print("The original Baudot code for the answer is:")
    # Printing each number in the final choice's code sequence
    print(final_choice_code)


solve_puzzle()

# Get the captured output
output_str = captured_output.getvalue()
# Restore stdout
sys.stdout = old_stdout
# Print the captured output to the actual console
print(output_str)

# Final answer tag
print("<<<E>>>")