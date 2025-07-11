# The final code is constructed by replacing @ with the 32-character solution.
# Note: The original template in the problem has a missing parenthesis. 
# The following is the corrected, working one-liner.
from re import *; print(*(sorted({int(n) for n in findall(r'\d+', input())}) or ["NO"]))