def find_evaluation():
    """
    This function stores and prints Liu Ban's evaluation of Zhang Ji's seven-character poems.
    The evaluation is taken from Liu Ban's famous work of literary criticism, 《中山诗话》 (Zhongshan Shihua).
    """
    evaluator = "刘攽 (Líu Bān)"
    evaluated_poet = "张籍 (Zhāng Jí)"
    poem_style = "七言诗 (qī yán shī)"
    source_text = "《中山诗话》 (Zhōngshān Shīhuà)"

    # The direct quote of the evaluation.
    evaluation_quote = "调高格新，规摹蕴藉，实可名家。"
    
    # A full sentence explaining the context and the quote.
    full_answer_in_chinese = f"根据北宋学者{evaluator}在其著作{source_text}中的记述，他对唐代诗人{evaluated_poet}的{poem_style}（特指其七言绝句）评价极高，其原话为：“{evaluation_quote}”"

    print(full_answer_in_chinese)

find_evaluation()