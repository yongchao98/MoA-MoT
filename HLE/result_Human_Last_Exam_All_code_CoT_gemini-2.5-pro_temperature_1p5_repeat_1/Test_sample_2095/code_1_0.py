def print_evaluation():
    """
    Prints Liu Ban's evaluation of Zhang Ji's seven-character poems.
    """
    source = "根据北宋刘攽在其诗歌评论著作《中山诗话》中的记载，他对张籍七言诗的评价如下。"
    original_quote_title = "原文为："
    original_quote = "“张籍七言诗，盖未为深妙，大概过于古质，非高流也。独其《行路难》一篇，稍优于此，可以比兴当时。”"
    
    explanation_title = "\n简而言之，刘攽的评价可以分解为以下几点："
    point_1 = "1. 整体评价不高：他认为张籍的七言诗整体上“未为深妙”，即不够深刻、精妙。"
    point_2 = "2. 主要批评：他认为其主要缺点在于风格“过于古质”，意思是诗歌的风格过于古朴、平实，缺乏文采。"
    point_3 = "3. 等级判定：他明确指出这些诗不属于上乘之作，即“非高流也”。"
    point_4 = "4. 特例赞扬：他特别指出了《行路难》这首诗是一个例外，认为它“稍优于此”（比其他的七言诗要好一些），并称赞其“可以比兴当时”，即能够用比兴的手法来讽喻或反映当时的社会现实。"

    print(source)
    print(original_quote_title)
    print(original_quote)
    print(explanation_title)
    print(point_1)
    print(point_2)
    print(point_3)
    print(point_4)

print_evaluation()