def get_authors():
    """
    Identifies and prints the authors of the 9 artworks.
    """
    authors = [
        "Liao Bingxiong",
        "Qi Baishi",
        "Jiang Zhaohe",
        "Dong Qichang",
        "Liu Guosong",
        "Zhao Mengfu",
        "Dong Qichang",
        "Yang Feiyun",
        "Ye Qianyu"
    ]

    print("The authors of the 9 works are, from left to right and top to bottom:")
    for i, author in enumerate(authors):
        print(f"Work {i+1}: {author}")

if __name__ == "__main__":
    get_authors()