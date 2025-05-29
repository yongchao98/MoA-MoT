def analyze_string(s):
    print(f"Analyzing string: {s}")
    print(f"Length: {len(s)}")
    print("\nCharacter counts:")
    for char in sorted(set(s)):
        print(f"'{char}': {s.count(char)}")
    
    print("\nChecking all rules:")
    print("1. 'aa' prefix?", s.startswith('aa'))
    print("2. Even number of 'b's?", s.count('b') > 0 and s.count('b') % 2 == 0)
    print("3. 'ab' prefix?", s.startswith('ab'))
    print("4. Contains 'acb'?", 'acb' in s)
    print("5. Contains 'ca' after first char?", 'ca' in s[1:])
    print("6. 'cb' prefix?", s.startswith('cb'))

s = "bbbcbbaaabbbb"
analyze_string(s)