# Our constructed code
print(*sorted({int(n) for n in findall(r'\d+', input())}) or ["NO"])

# Matching against the template's end part: `)})or ["NO"])`
# 1. The `)` closes `findall(...)`.
# 2. The `}` closes the set comprehension `{...}`.
# 3. The final `)` closes the `sorted(...)` call.
# The resulting `)})` sequence perfectly matches the template.