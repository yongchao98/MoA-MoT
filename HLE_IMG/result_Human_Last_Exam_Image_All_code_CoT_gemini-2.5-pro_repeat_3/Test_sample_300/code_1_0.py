def solve():
    """
    This function identifies and prints the authors of the 9 artworks.
    """
    authors = [
        "Huang Zhou (黄胄)",
        "Wu Guanzhong (吴冠中)",
        "Lin Sanzhi (林散之)",
        "Lin Sanzhi (林散之)",
        "Jin Shangyi (靳尚谊)",
        "Cheng Shifa (程十发)",
        "Wu Guanzhong (吴冠中)",
        "Li Keran (李可染)",
        "Qi Baishi (齐白石)"
    ]

    print("The authors of the 9 works, from left to right and top to bottom, are:")
    for i, author in enumerate(authors, 1):
        print(f"{i}. {author}")

solve()