def is_palindrome(s, left, right):
    while left < right:
        if s[left] != s[right]:
            return False
        left += 1
        right -= 1
    return True

def partition_palindrome(s):
    def backtrack(start, path):
        if start == len(s):
            result.append(path[:])
            return
        for end in range(start, len(s)):
            if is_palindrome(s, start, end):
                path.append(s[start:end+1])
                backtrack(end+1, path)
                path.pop()

    result = []
    backtrack(0, [])
    return result

s = "ddhahjhhjvyyv"
partitions = partition_palindrome(s)
print(partitions)