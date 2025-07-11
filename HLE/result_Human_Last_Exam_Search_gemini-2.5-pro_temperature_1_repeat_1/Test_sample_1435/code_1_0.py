def vigenere_subtract(s1, s2):
    res = ""
    for i in range(len(s1)):
        c1 = ord(s1[i]) - ord('a')
        c2 = ord(s2[i]) - ord('a')
        res += chr(((c1 - c2 + 26) % 26) + ord('a'))
    return res

def reverse_string(s):
    return s[::-1]

p = []
p.append("zuoeswzgnadou")  # P_1000
p.append("jqhonagacyrqj")  # P_999

# We need to find P_1, which is the 999th element before P_1000
# p_list[i] corresponds to P_{1000-i}
# We need p_list[999] = P_1

for i in range(2, 999 + 1):
    p_next = reverse_string(vigenere_subtract(p[i-2], p[i-1]))
    p.append(p_next)

# The answer is the last element computed.
# p[2] = P_998, ..., p[999] = P_1