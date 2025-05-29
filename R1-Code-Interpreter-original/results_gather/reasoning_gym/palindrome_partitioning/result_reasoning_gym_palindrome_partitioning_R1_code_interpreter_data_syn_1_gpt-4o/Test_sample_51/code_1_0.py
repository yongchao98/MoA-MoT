def is_palindrome(s):
    return s == s[::-1]

def partition(s):
    def backtrack(start, path):
        if start == len(s):
            result.append(path[:])
            return
        for end in range(start + 1, len(s) + 1):
            if is_palindrome(s[start:end]):
                path.append(s[start:end])
                backtrack(end, path)
                path.pop()

    result = []
    backtrack(0, [])
    return result

s = "cjjcueecbcdd"
partitions = partition(s)

# Limit the number of partitions printed to avoid truncation
max_display = 10
for i, p in enumerate(partitions):
    if i < max_display:
        print(p)
    else:
        print(f"...and {len(partitions) - max_display} more partitions")
        break