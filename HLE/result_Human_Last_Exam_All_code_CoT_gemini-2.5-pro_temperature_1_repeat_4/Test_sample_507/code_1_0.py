def find_similar_words_in_asian_languages():
    """
    Identifies and displays surprising word similarities between two
    geographically and culturally distant Asian languages: Korean and Tamil.
    """
    language1 = "Korean"
    language2 = "Tamil"

    print(f"The two languages are {language1} (from East Asia) and {language2} (from South Asia).")
    print("These languages are from completely different language families and their cultures had no direct contact.\n")
    print("Here are the comparisons for the words 'mom', 'dad', and 'broom' (with a necessary clarification).\n")
    print("-" * 25)

    # 1. Mom
    print("\nWord: Mom\n")
    korean_mom = "eomma (엄마)"
    tamil_mom = "ammā (அம்மா)"
    print("The final equation for 'Mom' is:")
    print(f"In {language1}: {korean_mom}")
    print(f"In {language2}:  {tamil_mom}")
    print("Result: The words are nearly identical.\n")
    print("-" * 25)

    # 2. Dad
    print("\nWord: Dad\n")
    korean_dad = "appa (아빠)"
    tamil_dad = "appā (அப்பா)"
    print("The final equation for 'Dad' is:")
    print(f"In {language1}: {korean_dad}")
    print(f"In {language2}:  {tamil_dad}")
    print("Result: The words are identical.\n")
    print("-" * 25)

    # 3. Broom (interpreted as Grass)
    print("\nWord: Broom\n")
    print("The words for a manufactured 'broom' do not match.")
    print("However, traditional brooms are often made of grass. The words for 'grass' are identical.\n")
    korean_grass = "pul (풀)"
    tamil_grass = "pul (புல்)"
    print("The final equation for the material of a 'Broom' (Grass) is:")
    print(f"In {language1}, 'grass' is: {korean_grass}")
    print(f"In {language2}, 'grass' is:  {tamil_grass}")
    print("Result: The words are identical, suggesting a link via 'grass broom'.\n")
    print("-" * 25)


if __name__ == "__main__":
    find_similar_words_in_asian_languages()