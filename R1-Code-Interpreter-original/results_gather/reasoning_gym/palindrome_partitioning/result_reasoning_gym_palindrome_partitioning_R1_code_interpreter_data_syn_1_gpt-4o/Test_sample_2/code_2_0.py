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

s = "phpmycyyiiynn"
partitions = partition(s)

# Print the number of partitions and a sample of them
print(f"Total partitions: {len(partitions)}")
print("Sample partitions:")
for i, partition in enumerate(partitions[:5]):  # Display only the first 5 partitions for brevity
    print(f"{i + 1}: {partition}")