def print_poem():
  """
  Prints the full text of the poem "Spring View" by Du Fu.
  """
  title = "春望"
  author_info = "[唐] 杜甫"
  
  poem_lines = [
    "国破山河在，城春草木深。",
    "感时花溅泪，恨别鸟惊心。",
    "烽火连三月，家书抵万金。",
    "白头搔更短，浑欲不胜簪。"
  ]
  
  print(f"这幅书法作品是唐代诗人杜甫的《{title}》。\n")
  print("作品全文如下：\n")
  print(f"    《{title}》")
  print(f"    {author_info}\n")
  for line in poem_lines:
    print(line)

print_poem()