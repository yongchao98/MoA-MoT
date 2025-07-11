def ancient_chinese_poem_explanation():
    """
    Explains the calligraphy work and provides its full context in Simplified Chinese.
    """
    title = "《易经·乾卦·文言》 (节选)"
    author = "（周）姬昌"
    poem_lines = [
        "同声相应，同气相求。",
        "水流湿，火就燥。",
        "云从龙，风从虎。",
        "圣人作而万物睹。"
    ]
    
    explanation = [
        "你朋友的这幅书法作品是用鸟虫篆书写的，内容是“云从龙，风从虎”。",
        "这是一句非常著名的成语，源自《易经》。完整的原文内容是："
    ]
    
    print("\n".join(explanation))
    print("\n" + "="*30)
    print(f"      {title}")
    # print(f"        —— {author}") # Author is traditional, not strictly needed for the context.
    print("\n".join(poem_lines))
    print("="*30 + "\n")
    
    work_content_explanation = [
        "作品原文解析:",
        "这幅作品包含了这句名言，具体文字布局如下：",
        "右列（从上至下）：云 从 龙",
        "中列（从上至下）：风 从 虎",
        "书法家可能为了构图，在左列增加了“田 鸟 啸”作为意境的补充。",
        "因此，作品的完整内容可以理解为：云从龙，风从虎，田鸟啸。"
    ]
    print("\n".join(work_content_explanation))


ancient_chinese_poem_explanation()