# Helper functions to convert between text and integer lists
def text_to_int(text):
  return [ord(c) - ord('a') for c in text]

def int_to_text(ints):
  return "".join([chr(i + ord('a')) for i in ints])

# Vigenere subtraction
def vigenere_subtract(c, p):
  return [(ci - pi + 26) % 26 for ci, pi in zip(c, p)]

# Initial given strings
p1000_text = "zuoeswzgnadou"
e1000_text = "ikfcuwfgaoked"

# Convert to integer lists
p1000 = text_to_int(p1000_text)
e1000 = text_to_int(e1000_text)

# --- Step 1: Find P_999 ---
# rev(P_999) = E_1000 - P_1000
rev_p999_int = vigenere_subtract(e1000, p1000)
# P_999 is the reverse of rev(P_999)
p999 = rev_p999_int[::-1]

# --- Step 2: Iteratively find P_998, P_997, ..., P_1 ---
# Initialize the backward sequence with P_1000 and P_999
pk = p1000
pk_minus_1 = p999

# Loop 998 times to go from P_1000 down to P_2
for _ in range(998):
  # rev(P_{k-2}) = P_k - P_{k-1}
  rev_pk_minus_2 = vigenere_subtract(pk, pk_minus_1)
  # P_{k-2} is the reverse of this result
  pk_minus_2 = rev_pk_minus_2[::-1]

  # Move to the next step back
  pk = pk_minus_1
  pk_minus_1 = pk_minus_2

# The final pk_minus_1 is P_1
p1_int = pk_minus_1
p1_text = int_to_text(p1_int)

# p1_text would be the answer