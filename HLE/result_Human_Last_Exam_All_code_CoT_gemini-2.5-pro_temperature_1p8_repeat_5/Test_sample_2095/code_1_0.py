def print_evaluation():
    """
    This function prints Liu Ban's exact evaluation of Zhang Ji's seven-character poems,
    based on his work "Zhongshan Shihua" (《中山诗话》).
    """

    # The source is Liu Ban's "Zhongshan Shihua".
    source = "刘攽（北宋）在其《中山诗话》中对张籍的七言诗作出了如下评价："
    
    # The original quote of the evaluation.
    quote = "“张籍七言诗，拙多于巧，然吾常爱其‘海人无家海里住，采珠役象为岁赋。’又‘新沙场头一摊土，李娟张妹有主无？’语虽拙，意却工也。’”"
    
    # A brief explanation of the quote's meaning.
    explanation = "\n这段话的核心意思是：刘攽认为张籍的七言诗在语言上朴实无华（拙）多于华丽取巧（巧），但他高度赞赏诗歌虽然语言朴拙，其内在的意境和情感却非常精妙深刻（工）。"
    
    full_answer = source + "\n" + quote + explanation
    
    print(full_answer)

print_evaluation()