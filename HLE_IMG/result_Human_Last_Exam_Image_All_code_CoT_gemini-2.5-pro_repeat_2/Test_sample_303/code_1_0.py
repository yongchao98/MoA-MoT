def explain_calligraphy():
    """
    Prints the analysis and content of the provided calligraphy work.
    """
    title = "这幅书法作品的内容解析"
    separator = "-" * 30

    # The identified characters and phrases
    characters = {
        "phrase1": {"hanzi": "闻飞龙", "pinyin": "wén fēi lóng"},
        "phrase2": {"hanzi": "关山弹", "pinyin": "guān shān tán"},
        "phrase3": {"hanzi": "思鸟鸣", "pinyin": "sī niǎo míng"}
    }

    # Explanation text
    explanation = [
        "这幅作品是用一种非常艺术化的古代字体——鸟虫篆写成的。",
        "其阅读顺序是遵循中国传统，从右至左，每列从上到下阅读。",
        "作品包含九个字，可分为三句诗：",
        "",
        f"  {characters['phrase1']['hanzi']} ({characters['phrase1']['pinyin']})",
        f"  {characters['phrase2']['hanzi']} ({characters['phrase2']['pinyin']})",
        f"  {characters['phrase3']['hanzi']} ({characters['phrase3']['pinyin']})",
        "",
        "这几句诗很可能并非出自某一首著名的古诗，而更像是一幅现代书法家创作的集句作品，旨在营造一种古典的意境。",
        "其诗意可以解读为：",
        "“闻飞龙”：仿佛听见巨龙在天际飞舞，比喻英雄豪迈的志向或威严的气象。",
        "“关山弹”：在雄伟的边关山脉之间，或弹奏琴瑟以抒怀，或拉开弓弦以习武，描绘了边塞的生活场景。",
        "“思鸟鸣”：心中思念着家乡鸟儿清脆的鸣叫声，寄托了对和平、安宁生活的向往。",
        "",
        "总的来说，这幅作品通过组合经典的意象，描绘了一位边塞戍客的所闻、所为、所思，意境开阔，雄浑与柔情并存，是一件富有感染力的艺术品。"
    ]

    # Print the full explanation
    print(title)
    print(separator)
    for line in explanation:
        print(line)

if __name__ == '__main__':
    explain_calligraphy()
