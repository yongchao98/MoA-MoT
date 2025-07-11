def find_authors():
    """
    Identifies and prints the authors of the 9 artworks from the image.
    The order is from left to right, top to bottom.
    """
    authors = [
        "1. Huang Zhou (黄胄)",
        "2. Wu Guanzhong (吴冠中)",
        "3. Lin Sanzhi (林散之)",
        "4. Yu Youren (于右任)",
        "5. Jin Shangyi (靳尚谊)",
        "6. Lin Fengmian (林风眠)",
        "7. Wu Guanzhong (吴冠中)",
        "8. Li Keran (李可染)",
        "9. Qi Baishi (齐白石)"
    ]

    print("The authors of the 9 works are:")
    for author in authors:
        print(author)

if __name__ == "__main__":
    find_authors()