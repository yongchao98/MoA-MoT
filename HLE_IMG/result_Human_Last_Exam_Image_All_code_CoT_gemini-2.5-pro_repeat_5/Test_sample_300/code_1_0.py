def solve():
    """
    This function identifies and lists the authors of the 9 artworks.
    """
    authors = [
        "Huang Zhou (黄胄)",
        "Wu Guanzhong (吴冠中)",
        "Mao Zedong (毛泽东)",
        "Sha Menghai (沙孟海)",
        "Jin Shangyi (靳尚谊)",
        "Lin Fengmian (林风眠)",
        "Wu Guanzhong (吴冠中)",
        "Li Keran (李可染)",
        "Qi Baishi (齐白石)"
    ]

    print("The authors of the 9 works, from left to right, top to bottom are:")
    for i, author in enumerate(authors, 1):
        print(f"{i}. {author}")

solve()