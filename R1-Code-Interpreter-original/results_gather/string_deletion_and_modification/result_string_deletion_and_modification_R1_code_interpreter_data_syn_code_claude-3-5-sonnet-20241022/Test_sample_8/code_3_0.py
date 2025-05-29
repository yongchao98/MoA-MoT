def verify_string(s):
    print(f"Analyzing string: {s}")
    print(f"Length: {len(s)}")
    print(f"Ends with 'aa'?: {s.endswith('aa')}")
    print(f"Ends with 'bb'?: {s.endswith('bb')}")
    print(f"Starts with 'ab'?: {s.startswith('ab')}")
    print(f"Starts with 'bc'?: {s.startswith('bc')}")
    print(f"Contains 'abc'?: {'abc' in s}")
    print(f"Contains 'bca'?: {'bca' in s}")
    
    # Print all 3-character substrings
    print("\nAll 3-character substrings:")
    for i in range(len(s)-2):
        print(f"Position {i}: {s[i:i+3]}")

verify_string("cbccabac")