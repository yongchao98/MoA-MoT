def solve():
    """
    Identifies and prints the authors of the 9 artworks.
    """
    authors = [
        "Zhou Sicong (周思聪)",
        "Wu Guanzhong (吴冠中)",
        "Mao Zedong (毛泽东)",
        "Shi Lu (石鲁)",
        "Jin Shangyi (靳尚谊)",
        "Lin Fengmian (林风眠)",
        "Wu Guanzhong (吴冠中)",
        "Li Keran (李可染)",
        "Qi Baishi (齐白石)"
    ]

    print("From left to right, top to bottom, the authors are:")
    for i, author in enumerate(authors):
        print(f"{i+1}. {author}")

solve()
<<<1. Zhou Sicong, 2. Wu Guanzhong, 3. Mao Zedong, 4. Shi Lu, 5. Jin Shangyi, 6. Lin Fengmian, 7. Wu Guanzhong, 8. Li Keran, 9. Qi Baishi>>>