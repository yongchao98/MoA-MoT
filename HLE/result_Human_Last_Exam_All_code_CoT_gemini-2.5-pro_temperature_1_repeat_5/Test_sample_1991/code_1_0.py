def generate_descriptive_phrase():
    """
    This function holds and prints a fourteen-syllable phrase.
    The phrase describes abnormal cell growth in the Luys' body (subthalamic nucleus)
    and matches the scansion and rhyme of a famous Disney song title.
    """
    # The phrase is constructed to have 14 syllables and rhyme with "Supercalifragilisticexpialidocious".
    # Syllable count: The(1) Sub(4)-tha(4)-lam(4)-ic(4) Nu(3)-cle(3)-us(3) Growth(1) Is(1) So(1) A(3)-tro(3)-cious(3) = 14
    phrase = "The Subthalamic Nucleus Growth Is So Atrocious"
    
    # Capitalize each word, as per the instructions.
    capitalized_phrase = ' '.join([word.capitalize() for word in phrase.split()])
    
    print(capitalized_phrase)

generate_descriptive_phrase()