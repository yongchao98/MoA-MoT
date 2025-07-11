def explain_calligraphy():
    """
    This function provides an explanation for the provided calligraphy image,
    which is in the style of Xu Bing's "Book from the Sky".
    """
    title = "关于这幅书法作品的解释"
    line1 = "您好！您朋友展示的这幅书法作品，并非一首能够按字词解读的传统古诗。"
    line2 = "这件作品的风格，与中国当代艺术家徐冰先生创作的著名艺术项目《天书》(Tianshu, or Book from the Sky) 如出一辙。"
    line3 = "《天书》的核心特点是，其中所有的“汉字”都是艺术家本人精心设计、刻印出来的“伪汉字”。"
    line4 = "这些字在结构和笔画上都模仿了真实的汉字（您图中的风格类似鸟虫篆），让人感觉似曾相识，但实际上它们没有任何读音和意义，所以任何人都无法读懂。"
    line5 = "因此，您朋友说它“相对著名”，这是一个非常巧妙的说法。它的“著名”在于它作为一件挑战文字和文化定义的当代艺术品的概念，而不是因为它是一首大家耳熟能详的诗。"
    line6 = "所以，这幅作品的“完整内容”，就是它“无内容”这个概念本身。您的朋友跟您开了一个有趣的艺术玩笑！"

    print(f"{title}")
    print("-" * 30)
    print(f"{line1}\n")
    print(f"{line2}\n")
    print(f"{line3}\n")
    print(f"{line4}\n")
    print(f"{line5}\n")
    print(f"{line6}")

if __name__ == '__main__':
    explain_calligraphy()