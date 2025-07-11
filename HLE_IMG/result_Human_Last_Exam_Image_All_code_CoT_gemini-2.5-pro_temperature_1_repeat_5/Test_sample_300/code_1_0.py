def identify_authors():
    """
    Identifies and prints the authors of the 9 artworks.
    """
    authors = [
        "Huang Zhou (黄胄)",
        "Wu Guanzhong (吴冠中)",
        "Mao Zedong (毛泽东)",
        "Zhu Da / Bada Shanren (朱耷 / 八大山人)",
        "Jin Shangyi (靳尚谊)",
        "Lin Fengmian (林风眠)",
        "Wu Guanzhong (吴冠中)",
        "Guan Shanyue (关山月)",
        "Qi Baishi (齐白石)"
    ]

    print("The authors of the artworks from left to right, top to bottom are:")
    for i, author in enumerate(authors):
        print(f"{i+1}. {author}")

identify_authors()