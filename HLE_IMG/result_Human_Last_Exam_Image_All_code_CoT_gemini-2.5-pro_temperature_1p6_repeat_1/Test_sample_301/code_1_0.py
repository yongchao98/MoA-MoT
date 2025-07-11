def find_authors():
    """
    Identifies and prints the authors of the 9 artworks.
    """
    authors = [
        "Han Meilin",
        "Wu Changshuo",
        "Zhou Sicong",
        "Wang Xizhi",
        "Han Meilin",
        "Zhao Mengfu",
        "Mi Fu",
        "Jin Shangyi",
        "Huang Zhou"
    ]
    
    print("From left to right, top to bottom, the authors are:")
    for i, author in enumerate(authors, 1):
        print(f"{i}. {author}")
        
    final_answer = ", ".join(authors)
    print(f"\n<<<From left to right, from top to bottom, the author of each of the 9 works are: {final_answer}>>>")

find_authors()