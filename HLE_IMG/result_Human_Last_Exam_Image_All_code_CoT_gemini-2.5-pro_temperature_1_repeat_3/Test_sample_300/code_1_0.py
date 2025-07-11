def solve():
    """
    This function identifies and prints the authors of the nine artworks.
    """
    authors = [
        "Huang Zhou (黄胄)",
        "Li Keran (李可染)",
        "Mao Zedong (毛泽东)",
        "Yu Youren (于右任)",
        "Jin Shangyi (靳尚谊)",
        "Lin Fengmian (林风眠)",
        "Wu Guanzhong (吴冠中)",
        "Fu Baoshi (傅抱石)",
        "Qi Baishi (齐白石)"
    ]

    print("From left to right, from top to bottom, the authors are:")
    for i, author in enumerate(authors, 1):
        print(f"{i}. {author}")

solve()