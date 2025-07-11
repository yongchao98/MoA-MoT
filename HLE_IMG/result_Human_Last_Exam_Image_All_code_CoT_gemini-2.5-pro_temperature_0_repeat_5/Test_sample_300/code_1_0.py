def identify_authors():
    """
    Identifies and prints the authors of the nine artworks.
    """
    authors = [
        "1. Huang Zhou (黄胄)",
        "2. Wu Guanzhong (吴冠中)",
        "3. Mao Zedong (毛泽东)",
        "4. Xu Beihong (徐悲鸿)",
        "5. Jin Shangyi (靳尚谊)",
        "6. Lin Fengmian (林风眠)",
        "7. Wu Guanzhong (吴冠中)",
        "8. Li Keran (李可染)",
        "9. Qi Baishi (齐白石)"
    ]

    print("From left to right, top to bottom, the authors of the 9 works are:")
    for author in authors:
        print(author)

if __name__ == "__main__":
    identify_authors()