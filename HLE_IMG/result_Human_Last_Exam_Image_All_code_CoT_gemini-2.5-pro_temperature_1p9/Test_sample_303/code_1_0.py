def print_poem():
    """
    This function prints the full text of the famous Chinese poem "Bird-Chirping Stream" by Wang Wei,
    which is the subject of the calligraphy artwork provided.
    """
    
    explanation = "这幅书法作品的内容是唐代诗人王维的著名诗作《鸟鸣涧》。完整的诗歌内容如下："
    
    title = "《鸟鸣涧》"
    author = "（唐）王维"
    
    line1 = "人闲桂花落，"
    line2 = "夜静春山空。"
    line3 = "月出惊山鸟，"
    line4 = "时鸣春涧中。"

    print(explanation)
    print()  # for a blank line
    print(f"    {title}")
    print(f"    {author}")
    print()  # for a blank line
    print(line1)
    print(line2)
    print(line3)
    print(line4)

# Execute the function to display the poem
print_poem()