def sort(a, n):
  if n > 1:
    # First recursive call
    sort(a, n - 1)

    # Comparison and potential swap
    if a[n - 1] < a[n - 2]:
      swap(a[n - 1], a[n - 2])

      # Second recursive call (conditional)
      sort(a, n - 1)