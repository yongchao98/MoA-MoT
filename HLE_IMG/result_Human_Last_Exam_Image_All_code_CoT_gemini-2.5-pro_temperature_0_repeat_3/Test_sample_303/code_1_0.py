def print_poem():
    """
    This function prints the full content of the poem "Spring View" by Du Fu.
    The characters in the image correspond to the first three lines of this poem.
    """
    title = "《春望》"
    author = "杜甫 (唐)"
    poem = [
        "国破山河在，城春草木深。",
        "感时花溅泪，恨别鸟惊心。",
        "烽火连三月，家书抵万金。",
        "白头搔更短，浑欲不胜簪。"
    ]

    print("这幅书法作品是唐代诗人杜甫的著名诗作《春望》。")
    print("完整的诗歌内容如下（简体中文）：")
    print("\n")
    print(f"    {title}")
    print(f"      {author}\n")
    for line in poem:
        print(line)

print_poem()