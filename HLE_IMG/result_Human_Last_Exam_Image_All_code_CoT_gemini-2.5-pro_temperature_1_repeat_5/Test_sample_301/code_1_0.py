def get_artwork_authors():
    """
    Identifies and prints the authors of the 9 artworks in the image.
    The order is from left to right, top to bottom.
    """
    authors = [
        "1. Lao Shu (老树)",
        "2. Qi Baishi (齐白石)",
        "3. Ai Xuan (艾轩)",
        "4. Wang Xizhi (王羲之)",
        "5. Wu Guanzhong (吴冠中)",
        "6. Zhao Mengfu (赵孟頫)",
        "7. Huaisu (怀素)",
        "8. Jin Shangyi (靳尚谊)",
        "9. Lin Fengmian (林风眠)"
    ]

    print("From left to right, top to bottom, the authors are:")
    for author in authors:
        print(author)

if __name__ == "__main__":
    get_artwork_authors()