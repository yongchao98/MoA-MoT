import sys
import io

# Set stdout to utf-8 to handle Arabic characters
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

# 1. The name of the author
author_name = "maimonides"
print(author_name)

# Leave exactly one line space
print()

# 2. A list of the syllables with primary word stress for the first 10 words
# The first 10 words are:
# 1. أنا (anā) -> 'a-NĀ -> nā (نا)
# 2. أبيّن (ubayyinu) -> u-BAY-yi-nu -> bay (بي)
# 3. لك (laka) -> LA-ka -> la (ل)
# 4. أيضًا (ayḍan) -> 'AY-ḍan -> 'ay (أي)
# 5. في (fī) -> FĪ -> fī (في)
# 6. هذا (hādhā) -> HĀ-dhā -> hā (ها)
# 7. الفصل (al-faṣli) -> 'al-FAṢ-li -> faṣ (فص)
# 8. ودلائل (wa-dalā'ila) -> wa-da-LĀ-'i-la -> lā (لا)
# 9. التوحيد (al-tawḥīdi) -> 'al-taw-ḤĪ-di -> ḥī (حي)
# 10. على (ʿalā) -> 'a-LĀ -> lā (لى)
stressed_syllables = "نا بي ل أي في ها فص لا حي لى"
print(stressed_syllables)