from collections import Counter

def letter_frequency_analysis(ciphertext):
    # Remove non-alphabetic characters and convert to uppercase
    filtered_text = ''.join(filter(str.isalpha, ciphertext)).upper()
    # Count the frequency of each letter
    frequency = Counter(filtered_text)
    return frequency

ciphertext = "UGRF RFCGP FCYBQ ZCLR BMUL MTCP RFCGP CJCARPGA RFGPRW QAGCLRGDGA KCL UCPC YZQMPZCB GL RPYLQACLBCLRYJ AYJASJYRGMLQ."
frequency = letter_frequency_analysis(ciphertext)
print(frequency)