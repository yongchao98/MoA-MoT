def is_palindrome(s, left, right):
    while left < right:
        if s[left] != s[right]:
            return False
        left += 1
        right -= 1
    return True

def partition_helper(s, start, path, result):
    if start == len(s):
        result.append(path[:])
        return
    
    for end in range(start, len(s)):
        if is_palindrome(s, start, end):
            path.append(s[start:end+1])
            partition_helper(s, end+1, path, result)
            path.pop()

def partition(s):
    result = []
    partition_helper(s, 0, [], result)
    return result

# Given string
s = "phpmycyyiiynn"
# Get all palindrome partitions
partitions = partition(s)
print(partitions)