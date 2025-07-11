import textwrap

def solve():
    """
    Analyzes statements about Chinese wedding traditions to find the incorrect one.
    """
    analysis = {
        'A': "Correct. During the Ming Dynasty, the wedding attire for commoners began to emulate that of the nobility. This included a red robe, often with auspicious embroidery, and the phoenix crown (凤冠) and xiapei (霞帔). The red veil (红盖头) was also a standard part of the ensemble. The 'dragon and phoenix robe' (龙凤褂) is more specifically associated with the Qing dynasty and modern times, but the symbolism and general attire described are fitting for the Ming period's formalization of Han wedding dress.",
        'B': "Correct. According to ancient texts on rituals like the 'Book of Rites' (礼记), which describes Zhou dynasty customs, it was part of the etiquette for a bride to cover her head or face upon leaving her home, symbolizing her transition and modesty. This was a precursor to the later red veil.",
        'C': "Incorrect. While it was indeed a custom for brides in the Song Dynasty to use a fan to shield their face, the type of fan was typically a rigid, round or paddle-shaped fan (团扇, tuánshàn). The folding fan (折扇, zhéshàn) was introduced to China around this time but only became popular and widely used during the later Ming Dynasty. Therefore, specifying a folding fan for this period is an anachronism.",
        'D': "Correct. The xiapei (霞帔), a long, scarf-like vestment, was a mark of status for noblewomen starting in the Song Dynasty and became a standard part of formal wedding attire in the Ming Dynasty. The phoenix crown (凤冠), decorated with intricate designs of phoenixes and other birds using techniques like filigree and kingfisher feathers, was the headdress for empresses and high-ranking women, and its use was permitted for commoner brides on their wedding day during the Ming.",
        'E': "Correct. During the Qing Dynasty, Han Chinese wedding customs were largely preserved. The quintessential bridal color was red, symbolizing joy and good fortune. A typical outfit consisted of a red jacket or blouse (袄, ǎo) over a long, wide skirt (裙, qún), both heavily embroidered. Red embroidered shoes were also a standard part of the attire."
    }

    print("Analyzing each statement about Chinese wedding history:")
    print("-" * 50)
    for option, text in analysis.items():
        print(f"Statement {option}: {textwrap.fill(text, width=80)}\n")

    incorrect_option = 'C'
    print("-" * 50)
    print(f"Conclusion: Statement {incorrect_option} is the incorrect one because brides in the Tang and Song dynasties used rigid, round fans, not folding fans, to cover their faces.")

    print("<<<C>>>")

solve()