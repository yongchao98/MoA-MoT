def sort(a, n):
  if n > 1:
    # 1. Recursive call on the first n-1 elements
    sort(a, n - 1)

    # 2. Constant time operations (comparison and potential swap)
    if a[n - 1] < a[n - 2]:
      swap(a[n-1], a[n - 2])

    # 3. Another recursive call on the first n-1 elements
    sort(a, n - 1)