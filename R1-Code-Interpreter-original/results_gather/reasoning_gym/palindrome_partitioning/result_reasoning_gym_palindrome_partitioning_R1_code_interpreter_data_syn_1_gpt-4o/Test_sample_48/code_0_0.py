def is_palindrome(s):
    return s == s[::-1]

def partition(s):
    result = []
    part = []
    
    def backtrack(start):
        if start >= len(s):
            result.append(part.copy())
            return
        for end in range(start, len(s)):
            if is_palindrome(s[start:end+1]):
                part.append(s[start:end+1])
                backtrack(end+1)
                part.pop()
    
    backtrack(0)
    return result

s = "ossojnnja"
partitions = partition(s)
print(partitions)