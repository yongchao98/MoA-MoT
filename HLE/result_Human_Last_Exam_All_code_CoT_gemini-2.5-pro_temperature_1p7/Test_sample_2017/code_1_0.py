import sys
import io

# Ensure unicode characters are printed correctly
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

# Passage 1 Analysis Results
W1 = "εἶμαι"
C1 = "εἰμί"
P1 = "KoineByzantineDemotic"

# Passage 2 Analysis Results
W2 = "ὄρνιθα"
C2 = "ὄρνιν"
P2 = "KoineByzantineDemotic"

# Combine the results into the final specified format
final_answer = f"{W1},{C1},{P1},{W2},{C2},{P2}"

# Print the final answer string
print(final_answer)

# Present the answer in the requested final format
print(f'<<<{final_answer}>>>')