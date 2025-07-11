def find_authors():
    """
    Identifies and lists the authors of the 9 artworks.
    """
    authors = [
        "Xu Beihong (徐悲鸿)",
        "Wu Guanzhong (吴冠中)",
        "Qi Gong (启功)",
        "Xu Beihong (徐悲鸿)",
        "Jin Shangyi (靳尚谊)",
        "Cheng Shifa (程十发)",
        "Wu Guanzhong (吴冠中)",
        "Li Keran (李可染)",
        "Qi Baishi (齐白石)"
    ]

    print("From left to right, top to bottom, the authors are:")
    for i, author in enumerate(authors):
        print(f"Work {i+1}: {author}")

if __name__ == "__main__":
    find_authors()