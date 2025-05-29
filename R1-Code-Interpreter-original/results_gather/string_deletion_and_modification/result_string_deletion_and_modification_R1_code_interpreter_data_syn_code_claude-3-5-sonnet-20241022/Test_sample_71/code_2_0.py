def check_rules(s):
    print(f"Checking rules for string: {s}")
    # Rule 1: If prefix is 'ab'
    print(f"Rule 1 (ab prefix): {s.startswith('ab')}")
    # Rule 2: If starts with 'bb'
    print(f"Rule 2 (bb prefix): {s.startswith('bb')}")
    # Rule 3: If contains 'abc'
    print(f"Rule 3 (contains abc): {'abc' in s}")
    # Rule 4: If contains 'bca'
    print(f"Rule 4 (contains bca): {'bca' in s}")
    # Rule 5: If contains 'acb'
    print(f"Rule 5 (contains acb): {'acb' in s}")
    # Rule 6: If ends with 'aa'
    print(f"Rule 6 (ends with aa): {s.endswith('aa')}")

final_string = "bacccacaba"
check_rules(final_string)