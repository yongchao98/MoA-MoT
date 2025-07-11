def solve():
    """
    This function lists the authors of the 9 artworks shown in the image,
    ordered from left to right, top to bottom.
    """
    authors = [
        "Huang Zhou (黄胄)",       # Top-Left
        "Wu Guanzhong (吴冠中)",    # Top-Middle
        "Mao Zedong (毛泽东)",      # Top-Right
        "Qi Gong (启功)",          # Middle-Left
        "Jin Shangyi (靳尚谊)",      # Middle-Center
        "Lin Fengmian (林风眠)",     # Middle-Right
        "Wu Guanzhong (吴冠中)",    # Bottom-Left
        "Li Keran (李可染)",       # Bottom-Middle
        "Qi Baishi (齐白石)"       # Bottom-Right
    ]

    print("The authors of the 9 works, from left to right and top to bottom, are:")
    for i, author in enumerate(authors, 1):
        print(f"{i}. {author}")

solve()