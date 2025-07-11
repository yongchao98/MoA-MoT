# First, you may need to install the textblob library:
# pip install textblob
# Then, you may need to download the necessary data for it:
# python -m textblob.download_corpora

from textblob import TextBlob

def find_part_of_speech():
    """
    Analyzes a sentence to find the part of speech for the first word.
    """
    sentence = "Fat people eat accumulates."
    
    # Create a TextBlob object
    blob = TextBlob(sentence)
    
    # Get the part-of-speech tags for the sentence
    # This returns a list of (word, tag) tuples
    pos_tags = blob.tags
    
    # Get the first word and its tag
    first_word, first_tag = pos_tags[0]
    
    # Define the answer choices
    choices = {
        'A': 'Verb',
        'B': 'Noun',
        'C': 'Adjective',
        'D': 'Adverb',
        'E': 'Conjunction'
    }
    
    # Determine the full part of speech name from the tag
    # 'JJ' is the standard tag for an adjective
    pos_full_name = "Unknown"
    if first_tag.startswith('JJ'):
        pos_full_name = "Adjective"
    elif first_tag.startswith('NN'):
        pos_full_name = "Noun"
    elif first_tag.startswith('VB'):
        pos_full_name = "Verb"
    elif first_tag.startswith('RB'):
        pos_full_name = "Adverb"
        
    print(f"Sentence: \"{sentence}\"")
    print(f"The first word is: '{first_word}'")
    print(f"The Part-of-Speech (POS) tag identified by the library is: '{first_tag}'")
    print(f"The tag '{first_tag}' stands for: {pos_full_name}")

    # Find the matching choice
    final_choice_letter = None
    for letter, name in choices.items():
        if name == pos_full_name:
            final_choice_letter = letter
            break

    print(f"\nFinal Analysis Equation:")
    print(f"Word ('{first_word}') + Role (modifies 'people') => Part of Speech ({pos_full_name}) => Answer Choice ({final_choice_letter})")


find_part_of_speech()