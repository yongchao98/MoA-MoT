def get_evaluation_of_zhang_ji():
    """
    本函数用于展示北宋刘攽对唐代诗人张籍七言诗的评价。
    信息来源于刘攽的《中山诗话》。
    """
    
    # 定义评价的各个要素
    evaluator = "刘攽 (Líu Bǎn)"
    evaluated_poet = "张籍 (Zhāng Jí)"
    poem_style = "七言诗 (Seven-character poems)"
    
    # 评价原文的核心组成部分
    part_1 = "拙处似太白"  # 在质朴之处像李白（字太白）
    part_2 = "巧处似储光羲"  # 在工巧之处像储光羲
    conclusion = "所谓欲学好人，只见其似，未得其真者也" # 这就是所谓的想学习高手，却只学到皮毛，没得到精髓
    
    # 组合成完整的引文
    full_quote = f"{evaluated_poet}的{poem_style}，{part_1}，{part_2}，{conclusion}。"
    
    # 打印格式化的结果，模拟一个“评价等式”
    print("对张籍七言诗的评价分析如下：")
    print("==============================================")
    print(f"评价者: {evaluator}")
    print(f"被评价者: {evaluated_poet}")
    print("评价要点:")
    print(f"  (+) 质朴风格 ≈ 李白 ({part_1})")
    print(f"  (+) 工巧风格 ≈ 储光羲 ({part_2})")
    print(f"  (=) 最终结论: {conclusion}")
    print("==============================================")
    print("\n完整原文:")
    print(f"“张籍七言诗，拙处似太白，巧处似储光羲，所谓欲学好人，只见其似，未得其真者也。”")
    print("\n这段评价的意思是：刘攽认为张籍的七言诗虽然努力模仿前代大师，其质朴的风格有李白的神韵，工巧的方面也看得到储光羲的影子，但终究只是学到了形似，未能真正掌握其精神内核，模仿的痕迹较重。")

# 执行函数并输出结果
get_evaluation_of_zhang_ji()
