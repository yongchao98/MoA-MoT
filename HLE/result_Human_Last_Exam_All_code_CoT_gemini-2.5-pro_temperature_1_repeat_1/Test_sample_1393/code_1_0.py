import sys
import io

# Capture original stdout to restore it later if needed, though not strictly necessary here.
original_stdout = sys.stdout
# Create a string buffer to hold the output
output_buffer = io.StringIO()
# Redirect stdout to the buffer
sys.stdout = output_buffer

def decode_morse(morse_code_str):
    """Decodes a morse code string to text."""
    MORSE_CODE_DICT = {
        '.-': 'A', '-...': 'B', '-.-.': 'C', '-..': 'D', '.': 'E', '..-.': 'F',
        '--.': 'G', '....': 'H', '..': 'I', '.---': 'J', '-.-': 'K', '.-..': 'L',
        '--': 'M', '-.': 'N', '---': 'O', '.--.': 'P', '--.-': 'Q', '.-.': 'R',
        '...': 'S', '-': 'T', '..-': 'U', '...-': 'V', '.--': 'W', '-..-': 'X',
        '-.--': 'Y', '--..': 'Z', '-----': '0', '.----': '1', '..---': '2',
        '...--': '3', '....-': '4', '.....': '5', '-....': '6', '--...': '7',
        '---..': '8', '----.': '9', '.-.-.-': '.', '--..--': ',', '..--..': '?',
        '-..-.': '/', '-.--.': '(', '-.--.-': ')', '.-..-.': '"'
    }
    words = morse_code_str.split(' / ')
    decoded_message = []
    for word in words:
        decoded_word = ''
        letters = word.split(' ')
        for letter in letters:
            if letter in MORSE_CODE_DICT:
                decoded_word += MORSE_CODE_DICT[letter]
        decoded_message.append(decoded_word)
    return ' '.join(decoded_message)

def decode_baudot(baudot_code_str):
    """Decodes a Baudot (ITA2) code string to text. Assumes LTRS shift."""
    BAUDOT_LTRS_DICT = {
        '00011': 'A', '11001': 'B', '01110': 'C', '01001': 'D', '00001': 'E',
        '01101': 'F', '11010': 'G', '10100': 'H', '00110': 'I', '01011': 'J',
        '01111': 'K', '10010': 'L', '11100': 'M', '01100': 'N', '11000': 'O',
        '10110': 'P', '10111': 'Q', '01010': 'R', '00101': 'S', '10000': 'T',
        '00111': 'U', '11110': 'V', '10011': 'W', '11101': 'X', '10101': 'Y',
        '10001': 'Z', '00100': ' '
    }
    decoded_text = ''
    codes = baudot_code_str.split(' ')
    for code in codes:
        if code in BAUDOT_LTRS_DICT:
            decoded_text += BAUDOT_LTRS_DICT[code]
    return decoded_text

# Step 1: Decode the Morse Code question
morse_question = "... . .-.. . -.-. - / - .... . / -.-. --- .-. .-. . -.-. - / --- .-. .. --. .. -. / --- ..-. / - .... . / ..-. --- .-.. .-.. --- .-- .. -. --. / ... . -. - . -. -.-. . / .-..-. - .... . / ..-. .- -. / -.-. .... .. . ..-. -.--. - .. .- -. -.--.- - .- .. -. / .- -. -.. / -... .- -. -.. .. - ... / .- .-. . / .-- --- .-. - .... .-.. . ... ... / - --- / -- . -. - .. --- -. .-.-.- / .-- .. - .... / --- -. . / ... .. -. --. .-.. . / ... .-- --- .-. -.. --..-- / .. / -.-. .- -. / -... .-.. --- -.-. -.- / - .-. --- --- .--. ... / .- / -- .. .-.. .-.. .. --- -. .-.-.- .-..-. / - .... . / -.-. .... --- -.-. .. . / --. .. ...- . -. / -... . .-.. .-.. --- .-- / .. ... / . -. -.-. .-. -.-- .--. - . -.. / ..- ... .. -. --. / -... .- ..- -.. --- - / -.-. --- -.. ."
decoded_question = decode_morse(morse_question)
print("--- Step 1: Decoded Question ---")
print(f"The decoded Morse code is: {decoded_question}\n")

# Step 2: Analyze the question and provide context
print("--- Step 2: Analysis and Context ---")
print("The quote is from the famous Chinese play 'Single Sword Meeting' (单刀会), featuring the heroic general Guan Yu from the Three Kingdoms period.")
print("The question asks for the 'origin' of this sentence. While the character is historical, this specific dramatic line originates from theatrical adaptations.")
print("The answer choices are different styles of Chinese opera, encoded in Baudot code.\n")

# Step 3: Decode the Baudot answer choices
print("--- Step 3: Decoding Answer Choices ---")
choices = {
    'A': '01111 00111 01100 00100 10111 00111',
    'B': '10100 00111 00011 01100 11010 00100 11100 00001 00110 11101 00110',
    'C': '10101 00111 00001 00100 01011 00111',
    'D': '10101 00111 00100 01011 00111',
    'E': '01110 10100 00111 00011 01100 00100 01011 00111'
}

for choice, code in choices.items():
    decoded_choice = decode_baudot(code)
    print(f"Choice {choice}: {code} -> {decoded_choice}")

print("\n--- Step 4: Final Conclusion ---")
print("The story of 'Single Sword Meeting' is a representative masterpiece of Kunqu Opera (昆曲, Kūnqǔ).")
print("Therefore, 'KUN QU' is the correct origin of the theatrical sentence.")
print("The correct choice is A.\n")

print("--- Final Answer Breakdown ---")
correct_choice_code = choices['A'].split(' ')
correct_choice_chars = list(decode_baudot(choices['A']))

print("The final answer is derived by decoding choice A:")
for code, char in zip(correct_choice_code, correct_choice_chars):
    print(f"{code} -> '{char}'")

# Restore original stdout
sys.stdout = original_stdout
# Get the content from the buffer
output_content = output_buffer.getvalue()
# Close the buffer
output_buffer.close()
# Print the captured output to the actual console
print(output_content)
print("<<<A>>>")