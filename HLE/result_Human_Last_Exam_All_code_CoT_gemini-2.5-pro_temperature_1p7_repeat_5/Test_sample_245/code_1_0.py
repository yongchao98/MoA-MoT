import re

def solve():
    """
    This function processes a Russian sentence to find and list multi-syllable words
    where the accent is not on the first syllable.
    """
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."
    
    # List of multi-syllable words from the text where the accent is not on the first syllable.
    # The stress is identified manually: Шашлы́к, горе́лым, прибежа́л, сосе́д, дошёл, его́.
    # Words with stress on the first syllable: за́пах, вско́ре, ви́дно, до́ма.
    # One-syllable words: и, до.
    target_words_set = {"шашлык", "горелым", "прибежал", "сосед", "дошёл", "его"}
    
    # Clean and split the text into words
    # This regex splits by spaces and punctuation, keeping only alphabetic sequences.
    words = re.findall(r'\b[а-яА-ЯёЁ]+\b', text)
    
    result = []
    # Iterate through the words in their order of appearance
    for word in words:
        # Check if the lowercase version of the word is in our target set
        if word.lower() in target_words_set:
            result.append(word)
            
    # Print the comma-separated list
    print(",".join(result))

solve()
<<<Шашлык,горелым,прибежал,сосед,дошёл,его>>>