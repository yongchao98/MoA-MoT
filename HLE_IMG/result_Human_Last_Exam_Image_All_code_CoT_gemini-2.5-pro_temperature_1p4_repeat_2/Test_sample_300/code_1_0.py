def solve():
    """
    This function identifies and lists the authors of the 9 artworks.
    """
    authors = [
        "Xu Beihong (徐悲鴻)",
        "Wu Guanzhong (吳冠中)",
        "Xu Beihong (徐悲鴻)",
        "Xu Beihong (徐悲鴻)",
        "Jin Shangyi (靳尚谊)",
        "Lin Fengmian (林風眠)",
        "Wu Guanzhong (吳冠中)",
        "Li Keran (李可染)",
        "Qi Baishi (齊白石)"
    ]

    print("The authors of the 9 works, from left to right and top to bottom, are:")
    for i, author in enumerate(authors, 1):
        print(f"Work {i}: {author}")

solve()