def analyze_string(s):
    print(f"Analyzing string: {s}")
    print(f"Length: {len(s)}")
    print("\nSearching for 'ca':")
    
    # Print all occurrences of 'ca'
    for i in range(len(s)-1):
        if s[i:i+2] == 'ca':
            print(f"Found 'ca' at position {i}")
    
    # Count all characters
    char_count = {}
    for c in s:
        char_count[c] = char_count.get(c, 0) + 1
    print("\nCharacter counts:", char_count)
    
    # Check specific substrings
    print("\nSpecific checks:")
    print("First character:", s[0])
    print("Last two characters:", s[-2:])
    print("Contains 'bca':", 'bca' in s)
    print("First two characters:", s[:2])
    
    # Find all overlapping occurrences of 'ca'
    all_ca = []
    pos = -1
    while True:
        pos = s.find('ca', pos + 1)
        if pos == -1:
            break
        all_ca.append(pos)
    print("\nAll 'ca' positions:", all_ca)

s = "aaccbcbaacccbbab"
analyze_string(s)