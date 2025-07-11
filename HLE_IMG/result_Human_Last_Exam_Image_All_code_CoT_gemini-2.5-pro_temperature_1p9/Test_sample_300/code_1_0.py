def solve():
    """
    This function identifies and lists the authors of the nine artworks.
    """
    authors = [
        "Liu Bangchu (刘邦楚)",   # 1. Top Left: Figure painting of Yi minority women.
        "Wu Guanzhong (吴冠中)",  # 2. Top Middle: Abstract ink landscape.
        "Fu Shan (傅山)",        # 3. Top Right: Cursive script calligraphy.
        "Bada Shanren (八大山人)", # 4. Middle Left: Cursive script calligraphy.
        "Jin Shangyi (靳尚谊)",  # 5. Middle Center: Realist oil portrait of a young woman painting.
        "Lin Fengmian (林风眠)",  # 6. Middle Right: Stylized painting of a woman in traditional dress.
        "Wu Guanzhong (吴冠中)",  # 7. Bottom Left: Colorful painting of a farmyard with gourds.
        "Li Keran (李可染)",     # 8. Bottom Middle: Landscape painting with mountains and red-sailed boats.
        "Qi Baishi (齐白石)"      # 9. Bottom Right: Ink painting of a lobster.
    ]

    print("From left to right, top to bottom, the authors are:")
    for i, author in enumerate(authors, 1):
        print(f"{i}. {author}")

solve()