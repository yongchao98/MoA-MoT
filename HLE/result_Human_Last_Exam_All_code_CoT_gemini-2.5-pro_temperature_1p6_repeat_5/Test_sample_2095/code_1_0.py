# This script prints the exact evaluation by Liu Ban (刘攽) on Zhang Ji's (张籍) seven-character poems.

def print_evaluation():
  """
  Prints Liu Ban's evaluation of Zhang Ji's poetry in Chinese.
  """
  evaluation_summary = "根据北宋刘攽在其著作《中山诗话》中的记载，他对唐代诗人张籍的七言诗的评价可以概括为“拙多于巧”。"
  
  original_quote = "原文是：“张籍七言诗，拙多于巧，然吾常爱其‘海燕飞来广陵郡，江帆到日客心惊’，以为佳句也。”"
  
  explanation = """
  这段话的意思是：刘攽认为张籍的七言诗整体风格质朴，甚至有些笨拙，胜过了其工整巧妙之处（拙多于巧）。
  然而，他并非完全否定，反而非常欣赏其中的佳句，并特别举出“海燕飞来广陵郡，江帆到日客心惊”这两句，认为它们是难得的佳作。
  因此，刘攽的评价是一种带有欣赏的批评，认为张籍的诗歌于质朴天成中见其功力。
  """
  
  print(evaluation_summary)
  print(original_quote)
  print(explanation)

print_evaluation()