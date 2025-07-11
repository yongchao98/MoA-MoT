def solve():
    """
    This function prints the full text of the famous Chinese poem
    from which the calligraphic work is derived.
    """
    title = "《凉州词·其一》"
    author = "[唐] 王之涣"
    poem_lines = [
        "黄河远上白云间，",
        "一片孤城万仞山。",
        "羌笛何须怨杨柳，",
        "春风不度玉门关。"
    ]

    print("这幅书法作品中的文字选自唐代诗人王之涣的著名诗歌《凉州词》。")
    print("完整的诗歌内容如下：\n")
    print(title)
    print(author)
    print("\n".join(poem_lines))

solve()