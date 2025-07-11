import collections

def decode_morse(morse_code_string):
    """Decodes a morse code string into English."""
    morse_code_dict = {
        '.-': 'A', '-...': 'B', '-.-.': 'C', '-..': 'D', '.': 'E', '..-.': 'F',
        '--.': 'G', '....': 'H', '..': 'I', '.---': 'J', '-.-': 'K', '.-..': 'L',
        '--': 'M', '-.': 'N', '---': 'O', '.--.': 'P', '--.-': 'Q', '.-.': 'R',
        '...': 'S', '-': 'T', '..-': 'U', '...-': 'V', '.--': 'W', '-..-': 'X',
        '-.--': 'Y', '--..': 'Z',
        '.----': '1', '..---': '2', '...--': '3', '....-': '4', '.....': '5',
        '-....': '6', '--...': '7', '---..': '8', '----.': '9', '-----': '0',
        '.-.-.-': '.', '--..--': ',', '..--..': '?', '.----.': "'", '-.-.--': '!',
        '-..-.': '/', '-.--.': '(', '-.--.-': ')', '.-...': '&', '---...': ':',
        '-.-.-.': ';', '-...-': '=', '.-.-.': '+', '-....-': '-', '..--.-': '_',

        '.-..-.': '"'
    }

    decoded_message = ''
    words = morse_code_string.strip().split(' / ')
    for word in words:
        letters = word.split(' ')
        for letter in letters:
            if letter in morse_code_dict:
                decoded_message += morse_code_dict[letter]
        decoded_message += ' '
    return decoded_message.strip()

def decode_baudot(baudot_code_string):
    """Decodes a Baudot (ITA2) code string into English."""
    # Standard ITA2 Letters Shift (LTRS) table
    ita2_dict = {
        '00011': 'A', '11001': 'B', '01110': 'C', '01001': 'D', '00001': 'E',
        '01101': 'F', '11010': 'G', '10100': 'H', '00110': 'I', '01011': 'J',
        '01111': 'K', '10010': 'L', '11100': 'M', '01100': 'N', '11000': 'O',
        '10110': 'P', '10111': 'Q', '01010': 'R', '00101': 'S', '10000': 'T',
        '00111': 'U', '11110': 'V', '10011': 'W', '11101': 'X', '10101': 'Y',
        '10001': 'Z', '00100': ' '
        # Note: FIGS shift characters are ignored for this problem.
    }

    decoded_message = ''
    chars = baudot_code_string.strip().split(' ')
    for char in chars:
        if char in ita2_dict:
            decoded_message += ita2_dict[char]
    return decoded_message

# Step 1: Decode the Morse Code question
encoded_question = "... . .-.. . -.-. - / - .... . / -.-. --- .-. .-. . -.-. - / --- .-. .. --. .. -. / --- ..-. / - .... . / ..-. --- .-.. .-.. --- .-- .. -. --. / ... . -. - . -. -.-. . / .-..-. - .... . / ..-. .- -. / -.-. .... .. . ..-. -.--. - .. .- -. -.--.- - .- .. -. / .- -. -.. / -... .- -. -.. .. - ... / .- .-. . / .-- --- .-. - .... .-.. . ... ... / - --- / -- . -. - .. --- -. .-.-.- / .-- .. - .... / --- -. . / ... .. -. --. .-.. . / ... .-- --- .-. -.. --..-- / .. / -.-. .- -. / -... .-.. --- -.-. -.- / - .-. --- --- .--. ... / .- / -- .. .-.. .-.. .. --- -. .-.-.- .-..-. / - .... . / -.-. .... --- -.-. .. . / --. .. ...- . -. / -... . .-.. .-.. --- .-- / .. ... / . -. -.-. .-. -.-- .--. - . -.. / ..- ... .. -. --. / -... .- ..- -.. --- - / -.-. --- -.. ."
decoded_question = decode_morse(encoded_question)
print("--- Step 1: Decoded Question ---")
print(f"The Morse code translates to:\n'{decoded_question}'")
print("\n(Note: There appear to be minor typos like 'CHOCIE' for 'CHOICE' and 'BELLOW' for 'BELOW' in the original message, but the meaning is clear.)")

# Step 2: Analyze the quote
print("\n--- Step 2: Quote Analysis ---")
print("The quote contains the line: 'WITH ONE SINGLE SWORD, I CAN BLOCK TROOPS A MILLION.'")
print("This is a famous story from the classic Chinese novel 'Romance of the Three Kingdoms', describing the general Zhang Fei at the Battle of Changban.")
print("Zhang Fei was a hero of the Shu Kingdom, which was located in modern-day Sichuan province, China.")

# Step 3: Decode the Baudot answer choices
encoded_choices = {
    'A': "01111 00111 01100 00100 10111 00111",
    'B': "10100 00111 00011 01100 11010 00100 11100 00001 00110 11101 00110",
    'C': "10101 00111 00001 00100 01011 00111",
    'D': "10101 00111 00100 01011 00111",
    'E': "01110 10100 00111 00011 01100 00100 01011 00111"
}

print("\n--- Step 3: Decoded Answer Choices ---")
print("The answer choices are different forms of Chinese opera, encrypted in Baudot code. Let's decode them:")
decoded_results = {}
for choice, code in encoded_choices.items():
    decoded_text = decode_baudot(code)
    decoded_results[choice] = decoded_text
    print(f"Choice {choice}: {code} -> {decoded_text}")

# Step 4: Conclusion
print("\n--- Step 4: Conclusion ---")
print("The quote is about Zhang Fei, a hero from Sichuan.")
print("The decoded answers are all types of Chinese opera.")
print("The most logical origin for a story about a famous Sichuan hero is Sichuan Opera.")
print(f"Sichuan Opera is known as 'Chuan Ju', which matches the decoded text for choice E: '{decoded_results['E']}'.")
print("\nTherefore, the correct choice is E.")
final_choice_code = encoded_choices['E'].split(' ')
final_choice_text = list(decoded_results['E'])

print("\nThe final equation is:")
# Remove space for the text list
final_choice_text_no_space = [c for c in final_choice_text if c != ' ']
equation_parts = []
for i in range(len(final_choice_text_no_space)):
    equation_parts.append(f"{final_choice_code[i]} ({final_choice_text_no_space[i]})")

# Handle the space character
space_index = final_choice_text.index(' ')
print(f"{' '.join(equation_parts[:space_index])} + {final_choice_code[space_index]} (space) + {' '.join(equation_parts[space_index:])} = CHUAN JU")