import textwrap

def analyze_wedding_customs():
    """
    Analyzes statements about wedding customs in different Chinese dynasties to find the incorrect one.
    """
    # Dictionary holding the analysis for each statement.
    analysis = {
        'A': "CORRECT. The Ming Dynasty did indeed formalize much of what is considered traditional wedding attire. The 'Phoenix Crown and Xiapei' (凤冠霞帔) became the standard formal dress for noblewomen, and its elements, along with the red veil, were adopted as the ideal for bridal wear.",
        'B': "CORRECT. Ancient ceremonial texts, such as the 'Book of Rites' from the Zhou Dynasty, describe the bride's head being covered during the wedding procession for modesty and ritual purposes.",
        'C': "INCORRECT. The custom of a bride covering her face (障面, zhàngmiàn) when leaving home was practiced during the Tang and Song dynasties. However, the item used was typically a round fan (团扇, tuánshàn). The folding fan (折扇, zhéshàn) was introduced to China later and was primarily associated with male scholars, not bridal customs of that era.",
        'D': "CORRECT. Following the Song dynasty, especially in the Ming and Qing, the 'xiapei' (a decorative shawl) and 'phoenix crown' (凤冠) became the established formal attire for high-ranking women, bestowed by imperial decree. This elite style was emulated in wedding ceremonies.",
        'E': "CORRECT. During the Qing Dynasty, Han Chinese brides were generally permitted to wear traditional Ming-style clothing for their weddings. The ensemble of a red robe over a blouse, a red skirt, and embroidered shoes is a classic and accurate description of this attire."
    }

    # Print the analysis for each option
    print("Analyzing each statement about Chinese wedding customs:\n")
    for option, text in analysis.items():
        # Use textwrap for clean printing
        wrapped_text = textwrap.fill(text, width=80)
        print(f"Statement {option}:\n{wrapped_text}\n")

    # State the final conclusion
    print("Conclusion: Statement C contains a factual error regarding the type of fan used.")
    print("\nThe incorrect statement is C.")

# Execute the analysis
analyze_wedding_customs()

# The following line provides the final answer in the required format.
print("\n<<<C>>>")
