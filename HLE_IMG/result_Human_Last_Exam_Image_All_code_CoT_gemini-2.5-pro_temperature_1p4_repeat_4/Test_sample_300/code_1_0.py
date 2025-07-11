def solve():
    """
    Identifies and prints the authors of the 9 artworks.
    """
    authors = [
        "Zhou Sicong (周思聪)",
        "Wu Guanzhong (吴冠中)",
        "Xu Beihong (徐悲鸿)",
        "Xu Beihong (徐悲鸿)",
        "Jin Shangyi (靳尚谊)",
        "Lin Fengmian (林风眠)",
        "Wu Guanzhong (吴冠中)",
        "Li Keran (李可染)",
        "Qi Baishi (齐白石)"
    ]

    print("The authors of the artworks are:")
    for i, author in enumerate(authors):
        print(f"Work {i+1}: {author}")

solve()