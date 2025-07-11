def identify_artists():
    """
    Identifies and prints the author for each of the 9 artworks.
    """
    authors = [
        "Li Kuizheng (李魁正)",
        "Wu Guanzhong (吴冠中)",
        "Lin Sanzhi (林散之)",
        "Qi Gong (启功)",
        "Jin Shangyi (靳尚谊)",
        "Lin Fengmian (林风眠)",
        "Wu Guanzhong (吴冠中)",
        "Li Keran (李可染)",
        "Qi Baishi (齐白石)"
    ]

    print("The authors of the 9 works are:")
    for i, author in enumerate(authors, 1):
        print(f"{i}. {author}")

identify_artists()