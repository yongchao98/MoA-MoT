def check_possible_rules(s):
    print(f"Checking string: {s}")
    print(f"1. Contains 'abc': {'abc' in s}")
    print(f"2. Contains 'bca': {'bca' in s}")
    print(f"3. Even b's: {s.count('b')} b's")
    print(f"4. Starts with 'cc': {s.startswith('cc')}")
    print(f"5. Starts with 'ca': {s.startswith('ca')}")
    print(f"6. Starts with 'bc': {s.startswith('bc')}")

s = "accabbbcbbbaaaab"
check_possible_rules(s)