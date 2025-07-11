import textwrap

def solve_history_question():
    """
    Analyzes historical statements about Chinese weddings and identifies the incorrect one.
    """
    title = "Analysis of Wedding Customs in Chinese Dynasties"
    print(title)
    print("=" * len(title))

    analysis = {
        'A': "The Ming Dynasty did formalize elaborate red wedding attire, including the phoenix crown and red veil. This is largely accurate.",
        'B': "Covering the bride's head with a veil was a custom in the ancient Zhou Dynasty, as documented in ritual texts. This is accurate.",
        'C': "Brides in the Tang and Song dynasties did use fans to cover their faces. However, they used rigid, non-folding fans (團扇). The *folding fan* (摺扇) only became popular later in the Ming Dynasty. Therefore, this statement is incorrect.",
        'D': "The xiapei (a shawl) and the phoenix crown were indeed standard, high-status ceremonial attire for married women after the Song Dynasty. This is accurate.",
        'E': "During the Qing Dynasty, Han Chinese brides were permitted to wear traditional Ming-style red wedding garments, as described. This is accurate."
    }

    print("\nEvaluating each option:")
    for option, text in analysis.items():
        print(f"Statement {option}:")
        print(textwrap.fill(text, width=80))
        print("-" * 20)

    incorrect_statement_explanation = f"""
Conclusion:
The incorrect statement is C because it specifies the use of a 'folding fan' during the Tang and Song dynasties.
The type of fan used for the face-covering custom during that period was a rigid, non-folding fan, not a folding fan.
"""
    print(incorrect_statement_explanation)

solve_history_question()
print("<<<C>>>")