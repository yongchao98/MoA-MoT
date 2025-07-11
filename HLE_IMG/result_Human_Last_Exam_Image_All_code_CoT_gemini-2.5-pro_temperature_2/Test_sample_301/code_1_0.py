def solve_art_authors():
    """
    Identifies and prints the author of each of the 9 artworks in the grid,
    ordered from left to right, top to bottom.
    """
    authors = [
        "1. Ding Cong (丁聪)",
        "2. Wu Changshuo (吴昌硕)",
        "3. Huang Zhou (黄胄)",
        "4. Mao Zedong (毛泽东)",
        "5. Wu Guanzhong (吴冠中)",
        "6. Qi Gong (启功)",
        "7. Huaisu (怀素)",
        "8. Jin Shangyi (靳尚谊)",
        "9. Lin Fengmian (林风眠)"
    ]

    print("From left to right, top to bottom, the authors of the 9 works are:")
    for author_info in authors:
        print(author_info)

solve_art_authors()