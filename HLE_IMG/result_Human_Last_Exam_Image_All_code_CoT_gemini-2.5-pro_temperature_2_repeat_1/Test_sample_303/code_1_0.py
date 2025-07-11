def display_full_poem():
    """
    This function prints the complete content of the poem "Bird-Chirping Stream" by Wang Wei.
    The calligraphy in the image is an excerpt from this poem.
    """
    title = "《鸟鸣涧》"
    author = "【唐】 王维"
    poem_lines = [
        "人闲桂花落，",
        "夜静春山空。",
        "月出惊山鸟，",
        "时鸣春涧中。"
    ]

    print("这幅书法作品节选自唐代诗人王维的名篇《鸟鸣涧》。")
    print("全诗内容如下：\n")
    print(f"    {title}")
    print(f"    {author}")
    print("--------------------")
    for line in poem_lines:
        print(line)
    print("--------------------")

display_full_poem()