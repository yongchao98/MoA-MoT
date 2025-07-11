def find_authors():
    """
    This function lists the identified authors of the 9 artworks.
    """
    authors = [
        "Ding Cong (丁聪)",
        "Wu Changshuo (吴昌硕)",
        "Ye Qianyu (叶浅予)",
        "Mi Fu (米芾)",
        "Liu Kuo-sung (刘国松)",
        "Chu Suiliang (褚遂良)",
        "Huaisu (怀素)",
        "Jin Shangyi (靳尚谊)",
        "Lin Fengmian (林风眠)"
    ]

    print("The authors of the 9 works, from left to right and top to bottom, are:")
    for i, author in enumerate(authors):
        print(f"{i + 1}. {author}")

if __name__ == "__main__":
    find_authors()