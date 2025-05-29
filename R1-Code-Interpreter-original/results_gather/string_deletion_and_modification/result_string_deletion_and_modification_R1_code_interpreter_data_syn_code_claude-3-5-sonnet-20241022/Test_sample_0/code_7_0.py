def analyze_string(s):
    print(f"\nAnalyzing string: {s}")
    print(f"Length: {len(s)}")
    print("Checking each character pair:")
    for i in range(len(s)-1):
        print(f"Position {i}: '{s[i:i+2]}'")
    
    # Specific checks
    print("\nRule conditions:")
    print(f"1. Ends with 'ac'? {s.endswith('ac')}")
    print(f"2. Contains 'bbcc'? {'bbcc' in s}")
    print(f"3. Contains 'abc'? {'abc' in s}")
    print(f"4. Starts with 'bc'? {s.startswith('bc')}")
    print(f"5. Contains 'ca' after pos 1?", end=" ")
    ca_pos = s.find('ca', 1)
    print(f"Yes at pos {ca_pos}" if ca_pos > 0 else "No")
    print(f"6. Ends with 'bb'? {s.endswith('bb')}")

# Analyze current string
s = "bbbbbbacbabbba"
analyze_string(s)