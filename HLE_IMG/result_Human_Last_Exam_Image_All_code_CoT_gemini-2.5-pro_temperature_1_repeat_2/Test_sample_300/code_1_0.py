def solve():
    """
    This function identifies and lists the authors of the nine artworks.
    """
    authors = [
        "Huang Zhou (黄胄)",
        "Wu Guanzhong (吴冠中)",
        "Mao Zedong (毛泽东)",
        "Xu Beihong (徐悲鸿)",
        "Jin Shangyi (靳尚谊)",
        "Lin Fengmian (林风眠)",
        "Wu Guanzhong (吴冠中)",
        "Li Keran (李可染)",
        "Qi Baishi (齐白石)"
    ]

    print("The authors of the artworks from left to right, top to bottom are:")
    for i, author in enumerate(authors, 1):
        print(f"Work {i}: {author}")

solve()
<<<Huang Zhou, Wu Guanzhong, Mao Zedong, Xu Beihong, Jin Shangyi, Lin Fengmian, Wu Guanzhong, Li Keran, Qi Baishi>>>