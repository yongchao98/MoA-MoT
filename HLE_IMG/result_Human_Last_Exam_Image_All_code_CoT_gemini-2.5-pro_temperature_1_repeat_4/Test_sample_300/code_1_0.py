def solve():
    """
    This function identifies and prints the author for each of the 9 artworks.
    """
    authors = [
        "Xu Beihong (徐悲鸿)",
        "Wu Guanzhong (吴冠中)",
        "Mao Zedong (毛泽东)",
        "Shitao (石涛)",
        "Jin Shangyi (靳尚谊)",
        "Lin Fengmian (林风眠)",
        "Wu Guanzhong (吴冠中)",
        "Li Keran (李可染)",
        "Qi Baishi (齐白石)"
    ]

    positions = [
        "1. Top-left",
        "2. Top-center",
        "3. Top-right",
        "4. Middle-left",
        "5. Middle-center",
        "6. Middle-right",
        "7. Bottom-left",
        "8. Bottom-center",
        "9. Bottom-right"
    ]

    print("The authors of the 9 artworks, from left to right and top to bottom, are:")
    for i in range(len(authors)):
        print(f"{positions[i]}: {authors[i]}")

solve()
<<<1. Xu Beihong, 2. Wu Guanzhong, 3. Mao Zedong, 4. Shitao, 5. Jin Shangyi, 6. Lin Fengmian, 7. Wu Guanzhong, 8. Li Keran, 9. Qi Baishi>>>