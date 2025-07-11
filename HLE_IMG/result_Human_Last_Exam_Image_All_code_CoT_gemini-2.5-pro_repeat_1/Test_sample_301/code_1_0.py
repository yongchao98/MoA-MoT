def solve():
    """
    This function identifies and prints the authors of the 9 artworks.
    """
    authors = [
        "Han Meilin (韩美林)",
        "Qi Baishi (齐白石)",
        "Luo Zhongli (罗中立)",
        "Mao Zedong (毛泽东)",
        "Wu Guanzhong (吴冠中)",
        "Chu Suiliang (褚遂良)",
        "Mao Zedong (毛泽东)",
        "Jin Shangyi (靳尚谊)",
        "Ding Shaoguang (丁绍光)"
    ]

    print("From left to right, top to bottom, the authors of the 9 works are:")
    for i, author in enumerate(authors):
        print(f"{i + 1}. {author}")

solve()