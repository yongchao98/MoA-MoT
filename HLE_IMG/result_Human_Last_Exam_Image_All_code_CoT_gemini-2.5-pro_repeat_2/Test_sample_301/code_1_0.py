def find_authors():
    """
    Identifies and prints the authors of the 9 artworks.
    """
    authors = [
        "Zhu Xuanxian",
        "Qi Baishi",
        "Luo Zhongli",
        "Mi Fu",
        "Wu Guanzhong",
        "Chu Suiliang",
        "Huaisu",
        "Jin Shangyi",
        "Lin Fengmian"
    ]
    
    print("The authors of the 9 works, from left to right and top to bottom, are:")
    for i, author in enumerate(authors, 1):
        print(f"{i}. {author}")

find_authors()