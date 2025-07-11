import sys
import io

# Ensure the output is encoded in UTF-8 to correctly display Chinese characters.
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

# The simplified Chinese characters identified from the seal.
# The phrase is 寿不听虑 (shòu bù tīng lǜ), which means "For a long life, do not listen to worries."
char1 = "寿"
char2 = "不"
char3 = "听"
char4 = "虑"

# Print the characters.
print(f"The simplified Chinese characters on the seal are: {char1}{char2}{char3}{char4}")