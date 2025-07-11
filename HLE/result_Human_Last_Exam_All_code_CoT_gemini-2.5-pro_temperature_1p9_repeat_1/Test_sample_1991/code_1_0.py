def generate_phrase():
    """
    Generates and prints a fourteen-syllable phrase about abnormal cell growth
    causing cancer, matching the scansion of a famous Disney song title.
    """
    # The phrase describes cancer in the subthalamic nucleus (Body of Luys).
    # "Subthalamic Mitosis" refers to cell division in that specific brain region.
    # "Growing Morbid And Grossest" describes the cancerous tumor formation.
    # The syllable count is 4 + 3 + 2 + 2 + 1 + 2 = 14.
    # It mimics the scansion and musicality of "Supercalifragilisticexpialidocious".
    phrase = "Subthalamic Mitosis Growing Morbid And Grossest"
    
    # Printing each word in the final equation as requested
    words = phrase.split()
    print(*words)

generate_phrase()