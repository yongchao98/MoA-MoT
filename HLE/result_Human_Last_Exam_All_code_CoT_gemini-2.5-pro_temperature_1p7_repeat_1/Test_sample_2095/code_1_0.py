def find_evaluation():
  """
  This function retrieves and prints Liu Ban's evaluation of Zhang Ji's poetry.
  """
  # The evaluator, from the Northern Song Dynasty.
  evaluator = "刘攽"

  # The poet being evaluated, from the Tang Dynasty.
  poet = "张籍"
  
  # The subject of the evaluation.
  subject = "七言诗"

  # Liu Ban's evaluation is found in his work "Zhongshan Shihua" (《中山诗话》).
  # The quote is "张籍七言诗，拙多于巧。"
  # This means Zhang Ji's seven-character poems possess more artlessness/sincerity (拙) than cleverness/artifice (巧).
  evaluation_text = "张籍七言诗，拙多于巧。"

  source = "《中山诗话》"

  print(f"根据{evaluator}在其著作{source}中的记载，他对{poet}{subject}的精准评价是：")
  print(f"“{evaluation_text}”")

find_evaluation()