def sort(a, n):
  if n > 1:
    # 1. Recursive call on a smaller problem
    sort(a, n - 1)

    # 2. Comparison and potential swap (constant time operations)
    if a[n - 1] < a[n - 2]:
      swap(a[n - 1], a[n - 2])

      # 3. Another recursive call if the swap occurred
      sort(a, n - 1)