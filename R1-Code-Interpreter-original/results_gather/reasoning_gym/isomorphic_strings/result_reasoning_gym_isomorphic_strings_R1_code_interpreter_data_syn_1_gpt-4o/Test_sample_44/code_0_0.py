def are_isomorphic(s1, s2):
    if len(s1) != len(s2):
        return False
    
    mapping_s1_to_s2 = {}
    mapping_s2_to_s1 = {}
    
    for char1, char2 in zip(s1, s2):
        if char1 in mapping_s1_to_s2:
            if mapping_s1_to_s2[char1] != char2:
                return False
        if char2 in mapping_s2_to_s1:
            if mapping_s2_to_s1[char2] != char1:
                return False
        
        mapping_s1_to_s2[char1] = char2
        mapping_s2_to_s1[char2] = char1
    
    return True

# Test the function with the given strings
result = are_isomorphic("nwswwz", "xdyddk")
print(result)