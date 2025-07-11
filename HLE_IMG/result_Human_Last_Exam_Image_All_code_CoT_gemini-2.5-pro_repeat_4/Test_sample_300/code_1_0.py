def solve():
    """
    Identifies and prints the authors of the 9 artworks.
    """
    authors = [
        "Huang Zhou (黄胄)",
        "Shi Lu (石鲁)",
        "Lin Sanzhi (林散之)",
        "Lin Sanzhi (林散之)",
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