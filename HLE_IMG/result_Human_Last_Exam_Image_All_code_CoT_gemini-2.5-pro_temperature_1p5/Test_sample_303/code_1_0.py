import sys
import io

# Set stdout to utf-8 to ensure Chinese characters are printed correctly.
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def get_calligraphy_content():
  """
  This function explains and returns the full content of the calligraphy piece.
  The work is not a linear poem but an artistic arrangement of four phrases.
  """
  # The four identified phrases in Simplified Chinese.
  phrase1 = "云龙"  # 雲龍 (yún lóng) - Cloud Dragon
  phrase2 = "雄风"  # 雄風 (xióng fēng) - Heroic Wind
  phrase3 = "明月"  # 明月 (míng yuè) - Bright Moon
  phrase4 = "飞出关" # 飛出關 (fēi chū guān) - To Fly Out of the Pass
  
  # The final explanation and content.
  full_content = (
      f"这幅书法作品的完整内容并非一首传统的线性诗歌，"
      f"而是一个由四个充满意象的词组艺术性地组合而成的集句作品。\n"
      f"这四个词组分别是：{phrase1}、{phrase2}、{phrase3}、{phrase4}。"
  )
  return full_content

print(get_calligraphy_content())