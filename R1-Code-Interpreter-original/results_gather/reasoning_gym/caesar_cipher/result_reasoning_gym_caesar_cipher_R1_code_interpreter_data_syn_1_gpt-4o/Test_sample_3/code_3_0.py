from collections import Counter

def find_most_common_letter(text):
    # Filter only alphabetic characters and convert to uppercase
    filtered_text = [char.upper() for char in text if char.isalpha()]
    # Count the frequency of each letter
    letter_counts = Counter(filtered_text)
    # Find the most common letter
    most_common_letter, _ = letter_counts.most_common(1)[0]
    return most_common_letter

cipher_text = "ESP NZAJCTRSE WLHD ZQ ESP AWLNP HSPCP JZF LCP WZNLEPO LWDZ RZGPCY HSLE JZF NLY OZ HTES ESTD HZCV."
most_common_letter = find_most_common_letter(cipher_text)

# Calculate the shift assuming the most common letter corresponds to 'E'
shift = (ord(most_common_letter) - ord('E')) % 26

# Decrypt the text using the calculated shift
decrypted_message = decrypt_caesar_cipher(cipher_text, shift)
print(decrypted_message)