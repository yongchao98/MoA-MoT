def is_palindrome(s):
    return s == s[::-1]

def backtrack(start, path, result, s):
    if start == len(s):
        result.append(path[:])
        return
    for end in range(start + 1, len(s) + 1):
        if is_palindrome(s[start:end]):
            path.append(s[start:end])
            backtrack(end, path, result, s)
            path.pop()

def partition(s):
    result = []
    backtrack(0, [], result, s)
    return result

# Partition the string "ddmmdq"
partitions = partition("ddmmdq")
print(partitions)