def get_evaluation():
    """
    This function stores and prints Liu Ban's evaluation of Zhang Ji's poetry.
    """
    source = "北宋的刘攽在其《中山诗话》中"
    evaluation_of = "对张籍的七言诗的精确评价是："
    quote = "“拙多于巧，然意地有可爱者。”"
    explanation = """
这句话的意思是：张籍的七言诗在风格和技巧上显得质朴笨拙多于华丽灵巧，但是其诗歌所蕴含的意境和立意却有非常值得欣赏和喜爱的地方。"""

    full_answer = f"""{source}，{evaluation_of}

{quote}
{explanation}"""
    print(full_answer)

get_evaluation()