def are_isomorphic(s1, s2):
    if len(s1) != len(s2):
        return False
    
    map_s1_to_s2 = {}
    map_s2_to_s1 = {}
    
    for char1, char2 in zip(s1, s2):
        if char1 in map_s1_to_s2:
            if map_s1_to_s2[char1] != char2:
                return False
        if char2 in map_s2_to_s1:
            if map_s2_to_s1[char2] != char1:
                return False
        
        map_s1_to_s2[char1] = char2
        map_s2_to_s1[char2] = char1
    
    return True

s1 = "ememihq"
s2 = "ththjua"
result = are_isomorphic(s1, s2)
print(result)