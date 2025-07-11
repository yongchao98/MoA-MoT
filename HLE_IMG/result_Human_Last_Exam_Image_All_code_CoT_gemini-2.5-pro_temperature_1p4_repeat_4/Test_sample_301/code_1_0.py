def find_authors():
    """
    This function lists the authors of the 9 artworks shown in the image.
    The order is from left to right, starting from the top row.
    """
    authors = [
        "Han Meilin (韩美林)",
        "Wu Changshuo (吴昌硕)",
        "Huang Zhou (黄胄)",
        "Zhao Mengfu (赵孟頫)",
        "Liu Guosong (刘国松)",
        "Li Yong (李邕)",
        "Mi Fu (米芾)",
        "Jin Shangyi (靳尚谊)",
        "Lin Yong (林墉)"
    ]

    print("From left to right, from top to bottom, the authors are:")
    for i, author in enumerate(authors):
        print(f"{i+1}. {author}")

if __name__ == "__main__":
    find_authors()