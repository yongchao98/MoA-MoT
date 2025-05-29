s = "ccbcaaabbcbcacaccc"

# Check all rules
print("Rule 1 ('acb' exists):", 'acb' in s)
print("Rule 2 (starts with 'aa'):", s.startswith('aa'))
print("Rule 3 (starts with 'ca'):", s.startswith('ca'))
print("Rule 4 (starts with 'bb'):", s.startswith('bb'))
print("Rule 5 (ends with 'bb'):", s.endswith('bb'))
print("Rule 6 (starts with 'bc'):", s.startswith('bc'))