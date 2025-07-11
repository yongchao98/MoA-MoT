def find_authors():
    """
    Identifies and prints the author for each of the 9 artworks,
    ordered from left to right, top to bottom.
    """
    authors = [
        "Shi Guoliang (史国良)",
        "Wu Guanzhong (吴冠中)",
        "Zhang Xu (张旭)",
        "Huang Binhong (黄宾虹)",
        "Jin Shangyi (靳尚谊)",
        "Lin Fengmian (林风眠)",
        "Wu Guanzhong (吴冠中)",
        "Li Keran (李可染)",
        "Qi Baishi (齐白石)"
    ]

    print("The authors of the 9 works are:")
    for i, author in enumerate(authors):
        # The grid is 3x3. We can calculate row and column.
        row = i // 3 + 1
        col = i % 3 + 1
        print(f"Work {i+1} (Row {row}, Col {col}): {author}")

find_authors()