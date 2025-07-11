def solve():
    """
    Identifies and lists the authors of the 9 artworks.
    """
    authors = [
        "Cui Zifan",      # 1. Top Left
        "Qi Baishi",      # 2. Top Middle
        "Huang Zhou",     # 3. Top Right
        "Huaisu",         # 4. Middle Left
        "Liu Kuo-sung",   # 5. Middle Middle
        "Wang Xizhi",     # 6. Middle Right
        "Mao Zedong",     # 7. Bottom Left
        "Chen Yifei",     # 8. Bottom Middle
        "Ding Shaoguang"  # 9. Bottom Right
    ]

    print("From left to right, top to bottom, the authors are:")
    for i, author in enumerate(authors):
        print(f"{i+1}. {author}")
    
    final_answer = ", ".join(authors)
    print(f"\n<<< {final_answer} >>>")

solve()