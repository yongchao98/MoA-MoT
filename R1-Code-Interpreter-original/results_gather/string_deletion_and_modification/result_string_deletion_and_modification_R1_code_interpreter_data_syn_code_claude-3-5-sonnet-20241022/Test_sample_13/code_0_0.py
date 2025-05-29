def apply_rules(s):
    # Returns (new_string, rule_applied)
    
    # Rule 1: Replace first 'abc' with 'cab'
    if 'abc' in s:
        pos = s.index('abc')
        return s[:pos] + 'cab' + s[pos+3:], 1
    
    # Rule 2: If ends with 'ca', remove last character
    if s.endswith('ca'):
        return s[:-1], 2
    
    # Rule 3: If starts with 'bb', remove second character
    if s.startswith('bb'):
        return s[0] + s[2:], 3
    
    # Rule 4: Replace first 'acb' with 'bca'
    if 'acb' in s:
        pos = s.index('acb')
        return s[:pos] + 'bca' + s[pos+3:], 4
    
    # Rule 5: If even number of 'b's (and at least one), append 'ab'
    b_count = s.count('b')
    if b_count > 0 and b_count % 2 == 0:
        return s + 'ab', 5
    
    # Rule 6: If ends with 'cc', replace with 'b' and prepend 'a'
    if s.endswith('cc'):
        return 'a' + s[:-2] + 'b', 6
    
    return s, 0

# Initial string
s = "aaacaabcbabccba"
iterations = []
iterations.append(s)

while True:
    new_s, rule = apply_rules(s)
    if rule == 0 or new_s == s:  # No rule applied or no change
        break
    s = new_s
    iterations.append(f"Rule {rule}: {s}")

print("Initial:", iterations[0])
for i in range(1, len(iterations)):
    print(iterations[i])
print("\nFinal string:", s)