import textwrap

def get_evaluation():
    """
    This function prints Liu Ban's (刘攽) evaluation of Zhang Ji's (张籍) seven-character poems.
    """
    evaluation_text = """
    根据北宋刘攽在其著作《中山诗话》中的记载，他对唐代诗人张籍的七言诗的评价是具体且有细微差别的。

    他的确切评价原文如下：
    “张籍长篇不如短篇，而七言尤不及五言。然‘秋风漠漠向昏黑’之句，备见唐人气象，虽名家亦不多得。”

    这段评价可以解读为两点：
    1.  总体上，刘攽认为张籍的七言诗不如他的五言诗写得好（“七言尤不及五言”）。
    2.  然而，他特别指出，尽管如此，张籍七言诗《秋思》中的“秋风漠漠向昏黑”这一句，却 superbly 展现了唐代诗歌的宏大气象，是连著名诗人都很难写出的佳句。

    因此，刘攽的评价是：张籍的七言诗整体成就不及其五言诗，但其中仍有名句足以流传，并体现了盛唐风貌。
    """
    print(textwrap.dedent(evaluation_text).strip())

get_evaluation()