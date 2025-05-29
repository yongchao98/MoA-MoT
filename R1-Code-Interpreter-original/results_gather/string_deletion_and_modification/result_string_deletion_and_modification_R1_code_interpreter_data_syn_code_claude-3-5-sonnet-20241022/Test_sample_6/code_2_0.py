s = "caaaccabbcccabaaaaa"

# Check all conditions
print("Rule 1 (bc prefix):", s.startswith('bc'))
print("Rule 2 (contains abc):", 'abc' in s)
print("Rule 3 (ac suffix):", s.endswith('ac'))
print("Rule 4 (ab prefix):", s.startswith('ab'))
print("Rule 5 (contains acb):", 'acb' in s)
print("Rule 6 (bb suffix):", s.endswith('bb'))