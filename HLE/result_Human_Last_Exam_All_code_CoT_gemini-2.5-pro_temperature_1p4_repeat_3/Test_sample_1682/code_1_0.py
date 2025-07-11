import textwrap

def solve_wedding_history_question():
    """
    Analyzes statements about Chinese wedding history to find the incorrect one.
    """

    explanation = """
    The incorrect statement is A. Here is a breakdown of why:

    The primary inaccuracy in statement A lies with the 'embroidered dragon and phoenix robe' (龙凤褂, lóng fèng guà). While the phoenix crown and red veil were indeed key elements of Ming Dynasty bridal wear, the 'dragon and phoenix robe' as we know it today is more a product of the Qing Dynasty and modern times, especially in Southern China.

    During the Ming Dynasty, the dragon was a powerful symbol strictly reserved for the emperor. It would have been a serious offense for a commoner's wedding attire to feature a dragon. A bride's formal robe would typically feature a phoenix, a symbol of the empress and high-ranking women, but not a dragon. The most formal wedding attire for a Ming Dynasty noblewoman consisted of the phoenix crown (凤冠) and a xiapei (霞帔), a long, decorative stole-like garment.

    The other statements are historically accurate:
    B. Veils have been used since ancient times (Zhou Dynasty) to signify a bride's modesty.
    C. Using a fan to shield the bride's face was a notable custom in the Tang and Song dynasties.
    D. The xiapei and phoenix crown became standardized formal attire for women of rank after the Song Dynasty, particularly in the Ming.
    E. The description of a red jacket-and-skirt ensemble is accurate for a Han Chinese bride in the Qing Dynasty.
    """

    # Printing the formatted explanation and the final answer.
    print(textwrap.dedent(explanation).strip())
    print("\nTherefore, the incorrect statement is:")
    print("<<<A>>>")

solve_wedding_history_question()