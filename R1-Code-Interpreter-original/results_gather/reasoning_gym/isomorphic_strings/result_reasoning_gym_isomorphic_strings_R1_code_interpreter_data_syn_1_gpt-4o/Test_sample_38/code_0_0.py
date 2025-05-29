def are_isomorphic(s1, s2):
    if len(s1) != len(s2):
        return False
    
    map_s1_to_s2 = {}
    map_s2_to_s1 = {}
    
    for char1, char2 in zip(s1, s2):
        if char1 in map_s1_to_s2:
            if map_s1_to_s2[char1] != char2:
                return False
        else:
            map_s1_to_s2[char1] = char2
        
        if char2 in map_s2_to_s1:
            if map_s2_to_s1[char2] != char1:
                return False
        else:
            map_s2_to_s1[char2] = char1
    
    return True

# Test the function with the given strings
s1 = "my"
s2 = "co"
result = are_isomorphic(s1, s2)
print(result)