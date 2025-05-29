import math

# Triangle ABC
a_abc, b_abc, c_abc = 5, 12, 13
s_abc = (a_abc + b_abc + c_abc) / 2
area_abc = math.sqrt(s_abc * (s_abc - a_abc) * (s_abc - b_abc) * (s_abc - c_abc))

# Triangle DEF
a_def, b_def, c_def = 8, 15, 17
s_def = (a_def + b_def + c_def) / 2
area_def = math.sqrt(s_def * (s_def - a_def) * (s_def - b_def) * (s_def - c_def))

# Ratio of areas
ratio = area_abc / area_def
print(ratio)